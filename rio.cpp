//
// Copyright 2010-2011,2014 Ettus Research LLC
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <uhd/types/tune_request.hpp>
//#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/thread.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <boost/thread/thread.hpp> 
#include <boost/chrono.hpp>
#include <iostream>
#include <fstream>
#include <csignal>
#include <complex>
#include <fftw3.h>
#include <time.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

namespace po = boost::program_options;

int get_secs_since_midnight()
{
  time_t s=time(NULL)%86400;
  return(s);
}

float hamming_window(int n, int k)
{
  return(0.54-0.46*cos(2.0*3.1415*((double)n)/((double)k)));
}

static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

template<typename samp_type> void recv_to_file(
    uhd::usrp::multi_usrp::sptr usrp,
    const std::string &cpu_format,
    const std::string &wire_format,
    size_t samps_per_buff,
    unsigned long long num_requested_samples,
    double time_requested = 0.0,
    bool bw_summary = false,
    bool stats = false,
    bool null = false,
    bool enable_size_map = false,
    bool continue_on_bad_packet = false,
    int n_freqs=NULL,
    double *freqs=NULL,
    int n_bins=NULL,
    int *freq_idx_table=NULL,
    double rate=NULL,
    double f0=NULL
){
    unsigned long long num_total_samps = 0;
    //create a receive streamer
    uhd::stream_args_t stream_args(cpu_format,wire_format);
    uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);

    uhd::rx_metadata_t md;
    std::vector<samp_type> buff(samps_per_buff);
    std::ofstream outfile;
    bool overflow_message = true;

    //setup streaming
    uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);

    stream_cmd.num_samps = size_t(num_requested_samples);
    stream_cmd.stream_now = true;
    stream_cmd.time_spec = uhd::time_spec_t();
    rx_stream->issue_stream_cmd(stream_cmd);

    boost::system_time start = boost::get_system_time();
    unsigned long long ticks_requested = (long)(time_requested * (double)boost::posix_time::time_duration::ticks_per_second());
    boost::posix_time::time_duration ticks_diff;
    boost::system_time last_update = start;
    unsigned long long last_update_samps = 0;

    typedef std::map<size_t,size_t> SizeMap;
    SizeMap mapSizes;

    fftwf_complex *in, *out;
    fftwf_plan p;
    int fftlen=buff.size();
    in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftlen);
    out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftlen);
    float *avg = (float *) malloc(sizeof(float) * fftlen);

    float *w = (float *) malloc(sizeof(float) * fftlen);

    for(int i=0; i<fftlen ; i++)
      w[i]=hamming_window(i,fftlen);

    int count=0;
    double amp_pwr=0.0;
    p = fftwf_plan_dft_1d(fftlen, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    time_t tnow;
    time_t tprev;
    time(&tprev);
    // unix seconds for start of day
    time_t day_s=tprev-(tprev%86400);
    char riofname[1024];
    char riosfname[1024];
    FILE **frio=(FILE **)malloc(sizeof(FILE *)*n_freqs);
    FILE *sfile;
    for(int j=0 ; j<n_freqs; j++){
      sprintf(riofname,"rio-%1.2f-%06d.txt",freqs[j]/1e6,day_s);
      frio[j]=fopen(riofname,"a");
    }
    
    // rio-30-today.txt

    sprintf(riosfname,"spec-%06d.bin",day_s);
    sfile=fopen(riosfname,"a");

    while(not stop_signal_called and (num_requested_samples != num_total_samps or num_requested_samples == 0)) {
        boost::system_time now = boost::get_system_time();

        size_t num_rx_samps = rx_stream->recv(in, buff.size(), md, 3.0, enable_size_map);


	for(int i=0; i<fftlen ; i++)
	{
	  in[i][0]=w[i]*in[i][0];
	  in[i][1]=w[i]*in[i][1];
	}
	fftwf_execute(p);	

	for(int i=0; i<fftlen ; i++)
	{
	  float pwra=in[i][0]*in[i][0] +in[i][1]*in[i][1];
	  float pwr=out[i][0]*out[i][0] +out[i][1]*out[i][1];
	  avg[i]+=pwr;
	  amp_pwr+=pwra;
	}
	count++;
	time(&tnow);
	if(tnow != tprev)
	{
	  float scale=float(count);
	  for(int i=0; i<fftlen ; i++)
	  {
	    avg[i]=avg[i]/scale;
	  }
	  double amp_rms=sqrt(amp_pwr/((double)fftlen*count));
	  for(int j=0 ; j < n_freqs ; j++)
	  {
	    double mean_pwr=0.0;
	    for(int i=0; i<n_bins ; i++)
	    {
	      mean_pwr+=avg[freq_idx_table[j*n_freqs+i]];
	    }
	    mean_pwr=mean_pwr/float(n_bins);

	    fprintf(frio[j],"%d  %1.4f\n",(int)tnow,10*log10(mean_pwr));
	    printf("Channel %d %1.2f MHz time %d count %d power %f (dB) rms voltage %1.2f \n",j,freqs[j]/1e6,(int)tnow,count,10*log10(mean_pwr),32768.0*amp_rms);
	    if(tnow % 10 == 0)
	    {
	      fflush(frio[j]);  // touch file at least every 10 s
	    }
	  }
	  // write fftlen, time, sample rate, and center frequency
	  fwrite(&fftlen,sizeof(uint64_t),1,sfile);
	  printf("%d\n",sizeof(uint64_t));
	  fwrite(&tprev,sizeof(uint64_t),1,sfile);
	  fwrite(&rate,sizeof(double),1,sfile);
	  fwrite(&f0,sizeof(double),1,sfile);
	  // then fft
	  fwrite(avg,sizeof(float),fftlen,sfile);
	  //	  if(tnow %  == 0){
	  //	    char fname[1024];
	    //	    sprintf(fname,"spec-%06d.bin",(int)tnow); // write spectrum
	    //	    FILE *fout=fopen(fname,"w");
	    
	    //	    fclose(fout);
	    //	  }

	  count=0;
	  for(int i=0; i<fftlen ; i++)
	  {
	    avg[i]=0.0;
	  }
	  amp_pwr=0.0;
	  time(&tprev);
	  // unix seconds for start of day
	  time_t s=tprev-(tprev%86400);
	  if(s!=day_s)
	  {
	    for(int j=0 ; j < n_freqs ; j++)
	    {
	      char cmd[1024];
	      fclose(frio[j]);
	      sprintf(cmd,"ln -s rio-%1.2f-%06d.txt rio-%1.2f-yesterday.txt",freqs[j]/1e6,day_s,freqs[j]/1e6);
	      system(cmd);


	      sprintf(riofname,"rio-%1.2f-%06d.txt",freqs[j]/1e6,s);
	      day_s=s;
	      printf("new day %d\n",day_s);
	      frio[j]=fopen(riofname,"a");

	      sprintf(cmd,"ln -s rio-%1.2f-%06d.txt rio-%1.2f-today.txt",freqs[j]/1e6,s,freqs[j]/1e6);
	      system(cmd);
	    }
	    fclose(sfile);
	    sprintf(riosfname,"spec-%06d.bin",day_s);
	    sfile=fopen(riosfname,"a");

	  }


	}

        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
            std::cout << boost::format("Timeout while streaming") << std::endl;
            break;
        }
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW){
            if (overflow_message) {
                overflow_message = false;
                std::cerr << boost::format(
                    "Got an overflow indication. Please consider the following:\n"
                    "  Your write medium must sustain a rate of %fMB/s.\n"
                    "  Dropped samples will not be written to the file.\n"
                    "  Please modify this example for your purposes.\n"
                    "  This message will not appear again.\n"
                ) % (usrp->get_rx_rate()*sizeof(samp_type)/1e6);
            }
            continue;
        }
        if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
            std::string error = str(boost::format("Receiver error: %s") % md.strerror());
            if (continue_on_bad_packet){
                std::cerr << error << std::endl;
                continue;
            }
            else
                throw std::runtime_error(error);
        }

        num_total_samps += num_rx_samps;

    }

    stream_cmd.stream_mode = uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
    rx_stream->issue_stream_cmd(stream_cmd);
    fftwf_destroy_plan(p);
    fftwf_free(in); 
    fftwf_free(out);
    for(int j=0;j<n_freqs;j++)
      fclose(frio[j]);

    fclose(sfile);
}

typedef boost::function<uhd::sensor_value_t (const std::string&)> get_sensor_fn_t;


int UHD_SAFE_MAIN(int argc, char *argv[]){
    uhd::set_thread_priority_safe();

    //variables to be set by po
    std::string args, file, type, ant, subdev, ref, wirefmt, riofreqs;
    size_t total_num_samps, spb;
    double rate, freq, gain, bw, total_time, setup_time, riobw;

    //setup the program options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "help message")
        ("args", po::value<std::string>(&args)->default_value(""), "multi uhd device address args")
        ("type", po::value<std::string>(&type)->default_value("float"), "sample type: double, float, or short")
        ("nsamps", po::value<size_t>(&total_num_samps)->default_value(0), "total number of samples to receive")
        ("duration", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("time", po::value<double>(&total_time), "(DEPRECATED) will go away soon! Use --duration instead")
        ("spb", po::value<size_t>(&spb)->default_value(1024), "samples per buffer")
        ("rate", po::value<double>(&rate)->default_value(25e6), "rate of incoming samples")
        ("riofreqs", po::value<std::string>(&riofreqs)->default_value("30e6,36e6,40e6,45e6"), "riometer frequencies")
        ("riobw", po::value<double>(&riobw)->default_value(200e3), "riometer bw")
        ("freq", po::value<double>(&freq)->default_value(35e6), "RF center frequency in Hz")
        ("gain", po::value<double>(&gain), "gain for the RF chain")
        ("ant", po::value<std::string>(&ant), "antenna selection")
        ("subdev", po::value<std::string>(&subdev)->default_value("A:A"), "subdevice specification")
        ("bw", po::value<double>(&bw), "analog frontend filter bandwidth in Hz")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "reference source (internal, external, mimo)")
        ("wirefmt", po::value<std::string>(&wirefmt)->default_value("sc16"), "wire format (sc8 or sc16)")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
        ("progress", "periodically display short-term bandwidth")
        ("stats", "show average bandwidth on exit")
        ("sizemap", "track packet size and display breakdown on exit")
        ("continue", "don't abort on a bad packet")
        ("skip-lo", "skip checking LO lock status")
        ("int-n", "tune USRP with integer-N tuning")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    //print the help message
    if (vm.count("help")) {
        std::cout << boost::format("UHD RX samples to file %s") % desc << std::endl;
        std::cout
            << std::endl
            << "This application streams data from a single channel of a USRP device to a file.\n"
            << std::endl;
        return ~0;
    }

    bool bw_summary = vm.count("progress") > 0;
    bool stats = vm.count("stats") > 0;
    bool null = vm.count("null") > 0;
    bool enable_size_map = vm.count("sizemap") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;


    int64_t fftlen=(int64_t)spb;

    double binw=((double)rate)/((double)fftlen);
    int n_bins=(int) (riobw/binw);

    double *freqs=(double *)malloc(sizeof(double)*fftlen);
    for(int i=0; i<fftlen/2 ; i++){
      freqs[i]=((double)i*rate)/fftlen + freq;
    }
    for(int i=0; i<fftlen/2 ; i++){
      freqs[i+spb/2]=(double)(-fftlen/2+(double)i)*rate/((double)fftlen) + freq;
    }
    for(int i=0; i<spb ; i++){
      printf("bin %d f %1.2f\n",i,freqs[i]/1e6);
    }
      printf("nbins %d\n",n_bins);
    char *token;
    char delim[3]=",";
    char tmp[1024];
    strcpy(tmp,(char *)riofreqs.c_str());
    token=strtok(tmp,delim);
    int fi=0;
    while(token != NULL)
    {
      printf("%d %s\n",fi,token);
      token=strtok(NULL,delim);
      fi++;
    }
    int n_freqs=fi;
    double *dfreqs=(double *)malloc(sizeof(double)*n_freqs);
    strcpy(tmp,(char *)riofreqs.c_str());
    token=strtok(tmp,delim);
    fi=0;
    while(token != NULL)
    {
      dfreqs[fi]=atof(token);
      printf("d %d %f\n",fi,dfreqs[fi]);
      token=strtok(NULL,delim);
      fi++;
    }

    int *freq_idx_table=(int *)malloc(sizeof(int)*n_freqs*n_bins);

    for(int i=0 ; i < n_freqs; i++)
    {
      int min_idx=0;
      double minv=100e6;
      int bini=0;
      for(int j=0 ; j<fftlen ; j++){
	if(abs(dfreqs[i]-freqs[j])<riobw)
	{
	  freq_idx_table[i*n_freqs+bini]=j;
	  bini=(bini+1)%n_bins;
	}
      }
      for(int j=0;j<n_bins;j++)
      {
	printf("freq %d %f idx %d\n",i,dfreqs[i],freq_idx_table[i*n_freqs+j]);
      }
    }

    if (enable_size_map)
        std::cout << "Packet size tracking enabled - will only recv one packet at a time!" << std::endl;

    //create a usrp device
    std::cout << std::endl;
    std::cout << boost::format("Creating the usrp device with: %s...") % args << std::endl;
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);

    //Lock mboard clocks
    usrp->set_clock_source(ref);

    //always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("subdev")) usrp->set_rx_subdev_spec(subdev);

    std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;

    //set the sample rate
    if (rate <= 0.0){
        std::cerr << "Please specify a valid sample rate" << std::endl;
        return ~0;
    }
    std::cout << boost::format("Setting RX Rate: %f Msps...") % (rate/1e6) << std::endl;
    usrp->set_rx_rate(rate);
    std::cout << boost::format("Actual RX Rate: %f Msps...") % (usrp->get_rx_rate()/1e6) << std::endl << std::endl;

    //set the center frequency
    if (vm.count("freq")) { //with default of 0.0 this will always be true
        std::cout << boost::format("Setting RX Freq: %f MHz...") % (freq/1e6) << std::endl;
        uhd::tune_request_t tune_request(freq);
        if(vm.count("int-n")) tune_request.args = uhd::device_addr_t("mode_n=integer");
        usrp->set_rx_freq(tune_request);
        std::cout << boost::format("Actual RX Freq: %f MHz...") % (usrp->get_rx_freq()/1e6) << std::endl << std::endl;
    }

    //set the rf gain
    if (vm.count("gain")) {
        std::cout << boost::format("Setting RX Gain: %f dB...") % gain << std::endl;
        usrp->set_rx_gain(gain);
        std::cout << boost::format("Actual RX Gain: %f dB...") % usrp->get_rx_gain() << std::endl << std::endl;
    }

    //set the IF filter bandwidth
    if (vm.count("bw")) {
        std::cout << boost::format("Setting RX Bandwidth: %f MHz...") % (bw/1e6) << std::endl;
        usrp->set_rx_bandwidth(bw);
        std::cout << boost::format("Actual RX Bandwidth: %f MHz...") % (usrp->get_rx_bandwidth()/1e6) << std::endl << std::endl;
    }

    //set the antenna
    if (vm.count("ant")) usrp->set_rx_antenna(ant);
    int test_time = 1;
    //    boost::this_thread::sleep(boost::posix_time::seconds(setup_time)); //allow for some setup time
    boost::this_thread::sleep_for(boost::chrono::seconds(1)); //allow for some setup time

    if (total_num_samps == 0)
    {
        std::signal(SIGINT, &sig_int_handler);
        std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
    }

#define recv_to_file_args(format) \
    (usrp, format, wirefmt, spb, total_num_samps, total_time, bw_summary, stats, null, enable_size_map, continue_on_bad_packet, n_freqs, dfreqs, n_bins, freq_idx_table, rate, freq)
    recv_to_file<std::complex<float> >recv_to_file_args("fc32");

    //finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

    return EXIT_SUCCESS;
}
