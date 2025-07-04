all:
	g++ -O3 -o rio rio.cpp -lm -lfftw3f -luhd -lboost_program_options -lboost_system -lboost_thread -lboost_chrono -lboost_timer
