#include "utilities.h"

#include <cassert>
#include <csignal>
#include <iostream>
#include <fstream>
#include <limits>
#include "ext/boost/lexical_cast.hpp"
#include <iomanip>
#include "../timer.h"

using namespace std;


#if OPERATING_SYSTEM == LINUX
static void exit_handler(int exit_code, void *hint);
#elif OPERATING_SYSTEM == OSX
static void exit_handler();
#include <mach/mach.h>
#elif OPERATING_SYSTEM == WINDOWS || OPERATING_SYSTEM == CYGWIN
#include <windows.h>
#include <psapi.h>
#endif

// See issue 469 for the reasons we chose this limit.
static char *memory_padding = new char[32 * 1024];

static void out_of_memory_handler();
static void signal_handler(int signal_number);


void register_event_handlers() {
    // When running out of memory, release some emergency memory and
    // terminate.
    set_new_handler(out_of_memory_handler);

    // On exit or when receiving certain signals such as SIGINT (Ctrl-C),
    // print the peak memory usage.
#if OPERATING_SYSTEM == LINUX
    on_exit(exit_handler, 0);
#elif OPERATING_SYSTEM == OSX
    atexit(exit_handler);
#elif OPERATING_SYSTEM == CYGWIN || OPERATING_SYSTEM == WINDOWS
    // nothing
#endif
    signal(SIGABRT, signal_handler);
    signal(SIGTERM, signal_handler);
    signal(SIGSEGV, signal_handler);
    signal(SIGINT, signal_handler);
#if OPERATING_SYSTEM != WINDOWS
    // This causes problems, see issue479.
    //signal(SIGXCPU, signal_handler);
#endif
}

#if OPERATING_SYSTEM == LINUX || OPERATING_SYSTEM == OSX
#if OPERATING_SYSTEM == LINUX
void exit_handler(int, void *) {
#elif OPERATING_SYSTEM == OSX
void exit_handler() {
#endif
    print_peak_memory(false);
}
#endif

void exit_with(ExitCode exitcode) {
    switch (exitcode) {
    case EXIT_PLAN_FOUND:
        cout << "Solution found." << endl;
        break;
    case EXIT_CRITICAL_ERROR:
        cerr << "Unexplained error occurred." << endl;
        break;
    case EXIT_INPUT_ERROR:
        cerr << "Usage error occurred." << endl;
        break;
    case EXIT_UNSUPPORTED:
        cerr << "Tried to use unsupported feature." << endl;
        break;
    case EXIT_UNSOLVABLE:
        cout << "Task is provably unsolvable." << endl;
        break;
    case EXIT_UNSOLVED_INCOMPLETE:
        cout << "Search stopped without finding a solution." << endl;
        break;
    case EXIT_OUT_OF_MEMORY:
        cout << "Memory limit has been reached." << endl;
        break;
    default:
        cerr << "Exitcode: " << exitcode << endl;
        ABORT("Unknown exitcode.");
    }
    exit(exitcode);
}

static void out_of_memory_handler() {
    assert(memory_padding);
    delete[] memory_padding;
    memory_padding = 0;
    cout << "Failed to allocate memory. Released memory buffer." << endl;
    exit_with(EXIT_OUT_OF_MEMORY);
}

void signal_handler(int signal_number) {
    // See glibc manual: "Handlers That Terminate the Process"
    static volatile sig_atomic_t handler_in_progress = 0;
    if (handler_in_progress)
        raise(signal_number);
    handler_in_progress = 1;
    print_peak_memory(false);
    cout << "caught signal " << signal_number << " -- exiting" << endl;
    signal(signal_number, SIG_DFL);
    raise(signal_number);
}

int get_peak_memory_in_kb(bool use_buffered_input) {
    // On error, produces a warning on cerr and returns -1.
    int memory_in_kb = -1;

#if OPERATING_SYSTEM == OSX
    unused_parameter(use_buffered_input);
    // Based on http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
    task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

    if (task_info(mach_task_self(), TASK_BASIC_INFO,
                  reinterpret_cast<task_info_t>(&t_info),
                  &t_info_count) == KERN_SUCCESS)
        memory_in_kb = t_info.virtual_size / 1024;
#elif OPERATING_SYSTEM == WINDOWS || OPERATING_SYSTEM == CYGWIN
    unused_parameter(use_buffered_input);
    // The file /proc/self/status is present under Cygwin, but contains no peak memory info.
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PROCESS_MEMORY_COUNTERS *>(&pmc), sizeof(pmc));
    memory_in_kb = pmc.PeakPagefileUsage / 1024;
#else
    ifstream procfile;
    if (!use_buffered_input) {
        procfile.rdbuf()->pubsetbuf(0, 0);
    }
    procfile.open("/proc/self/status");
    string word;
    while (procfile.good()) {
        procfile >> word;
        if (word == "VmPeak:") {
            procfile >> memory_in_kb;
            break;
        }
        // Skip to end of line.
        procfile.ignore(numeric_limits<streamsize>::max(), '\n');
    }
    if (procfile.fail())
        memory_in_kb = -1;
#endif

    if (memory_in_kb == -1)
        cerr << "warning: could not determine peak memory" << endl;
    return memory_in_kb;
}

void print_peak_memory(bool use_buffered_input) {
    cout << "Peak memory: " << get_peak_memory_in_kb(use_buffered_input) << " KB" << endl;
}


bool is_product_within_limit(int factor1, int factor2, int limit) {
    assert(factor1 >= 0 && factor1 <= limit);
    assert(factor2 >= 0 && factor2 <= limit);
    return factor2 == 0 || factor1 <= limit / factor2;
}

bool get_GA_patterns_from_file(std::vector<std::vector<int> > &all_pattern_col,bool disjoint,double mutation_rate,int pdb_max_size){
  if(stored_GA_patterns.size()==0){
    cout<<"No patterns stored,calling load_GA_Patterns_from_file"<<endl;
    load_GA_Patterns_from_file();
  }
  all_pattern_col.clear();//just in case this was previously populated
  all_pattern_col.resize(1);//just in case this was previously populated
    std::string line;
    std::string temp;
    unsigned found2;
    unsigned next_pattern_pos;
    bool found_PDB=false;
    //ifstream in(log_file.c_str());
    //cout<<"log_file:"<<log_file<<",g_plan_filename:"<<g_plan_filename;
    string problem_name_mod=g_plan_filename;
    problem_name_mod+=":";
    cout<<",problem_name_mod:"<<problem_name_mod<<endl;
    std::string disjoint_pattern("disjoint_patterns:,");
    if(disjoint){
      disjoint_pattern+="1";
    }
    else{
      disjoint_pattern+="0";
    }
    
    cout<<disjoint_pattern<<endl;
    std::string mutation_rate_string("mp:,");
    //mutation_rate_string+=boost::lexical_cast<std::string>(mutation_rate);
    //mutation_rate_string+=std::to_string(mutation_rate);
    ///std::ostringstream strs;
    std::ostringstream strs;strs << std::fixed << std::setprecision(7);strs<<mutation_rate;
    //strs << mutation_rate;
    //std::string str = strs.str();
    mutation_rate_string+=strs.str();
    mutation_rate_string+=",";
    cout<<"mutation_rate_string:"<<mutation_rate_string<<endl;
    
    std::string pdb_max_size_string("size:,");
    std::ostringstream strs2;strs2 << std::fixed;strs2<<pdb_max_size;
    pdb_max_size_string+=strs2.str();
    pdb_max_size_string+=",";
    cout<<"pdb_max_size_string:"<<pdb_max_size_string<<endl;

    for(size_t pattern=0;pattern<stored_GA_patterns.size();pattern++){
     line=stored_GA_patterns.at(pattern);
     if( line.find(problem_name_mod)!=string::npos&&line.find(disjoint_pattern)!=string::npos
  &&line.find(mutation_rate_string)!=string::npos&&line.find(pdb_max_size_string)!=string::npos){
       cout<<"line:"<<line<<endl;
       found_PDB=true;
  
       unsigned current_pos=line.find("]");
       int num_databases = std::count(line.begin(), line.end(), ']') - 1;//the first ] is for the heuristic number
       //cout<<"num_databases:"<<num_databases<<endl;
       all_pattern_col.resize(num_databases);


       current_pos=line.find("[",current_pos+1);
       //unsigned found2 = line.find("[",found+1);
       next_pattern_pos=line.find("-",current_pos+1);
       for(int i=0;i<num_databases;i++){
  //cout<<"reading database"<<i<<endl;
  while(next_pattern_pos>line.find_first_of(",]",current_pos+1)){
  //while(true)
    current_pos=line.find_first_of("0123456789",current_pos);//so it points to the next variable
    if(current_pos>next_pattern_pos){
      //cout<<"skipping empty database"<<endl;
      while(current_pos>next_pattern_pos){
        next_pattern_pos=line.find("-",next_pattern_pos+1);
        //cout<<"skipped database "<<i<<" because it is empty"<<endl;
        //cout<<all_pattern_col.at(i);cout<<endl;
        i++;
      }
      i--;//need to decrease i by one or it will be one ahead by the for statement
      break;
    }
    //cout<<"\tcurrent_pos:"<<current_pos;
    found2 = line.find_first_of(",]",current_pos);
    //cout<<",found2:"<<found2;
    if(line.find_first_of(",]",current_pos)==string::npos){
      //cout<<",finished with pattern:"<<i<<",string finished"<<endl;
      break;
    }
    temp=line.substr(current_pos,found2-current_pos);
    //cout<<",next var:"<<temp<<",";
    int temp2=boost::lexical_cast<int>(temp);
    all_pattern_col.at(i).push_back(temp2);
    // cout<<",last int added:"<<all_pattern_col.at(i).back();
    current_pos=found2;//current_pos not pointing to next , or ]
    if(line.find_first_of(",",found2)>next_pattern_pos){
      //cout<<",finished with pattern:"<<i<<",pattern finished"<<endl;
      current_pos=line.find("[",next_pattern_pos);
      next_pattern_pos=line.find("-",next_pattern_pos+1);
      //cout<<"next_pattern_pos:"<<next_pattern_pos<<endl;
      break;
    }
  }
  //cout<<"database:"<<i<<",read pattern:";
  //cout<<all_pattern_col.at(i);cout<<endl;
       } 
       //Now add to our time calculations how long did iPDB originally take to generate this PDBs
       current_pos=line.find("time:");
       current_pos=line.find_first_of("0123456789",current_pos);
       temp=line.substr(current_pos);
       //double temp2=boost::lexical_cast<double>(temp);
       //overall_original_pdbs_time+=temp2;
     }
   }
    if(!found_PDB){
      //cout<<"No existing gaPDB to read for current problem!, will return dummy heuristic!:"<<endl;
      //exit(o);
    }
    /*  else{
      cout<<"Original GAPDB time:"<<temp2<<mutation_rate_string<<disjoint_pattern<<",Original_pdbs_time:"<<overall_original_pdbs_time<<endl;
    }*/
    return found_PDB;
}

void load_GA_Patterns_from_file(){
  std::string line;

  string task2 = problem_name2;
  
  size_t found = task2.find(".");

  string problem_name_mod = task2.substr(0, found);
  problem_name_mod += ".dat";
  problem_name_mod = "/" + problem_name_mod;
  problem_name_mod = domain_name + problem_name_mod;
  problem_name_mod = "dat/" + problem_name_mod;
  cout<<"problem_name_mod = "<<problem_name_mod<<endl;

  ifstream in(problem_name_mod.c_str());
  Timer load_GA_from_file_timer;
  
  cout<<"Calling load_GA_Patterns_from_file"<<endl;
  cout<<"log_file:"<<problem_name_mod<<",g_plan_filename:"<<g_plan_filename<<endl;
  bool problem_found=false;
    if( in.is_open())
    {
      cout<<"is_open true"<<endl;
      while( getline(in,line) ){
	  //cout<<"line:"<<line<<endl;
	if( line.find(g_plan_filename)!=string::npos) {
	    //cout<<"inside the line"<<endl;
	    stored_GA_patterns.push_back(line);
	    problem_found=true;
	}
      }
      if(problem_found){
	cout<<"problem_found among stored GAs:"<<g_plan_filename<<endl;
      }
    }
    in.close();
    cout<<"stored_GA_patterns.size:"<<stored_GA_patterns.size()<<",time:"<<load_GA_from_file_timer()<<endl;
}
