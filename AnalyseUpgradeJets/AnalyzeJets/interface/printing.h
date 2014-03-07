#ifndef PRINTING
#define PRINTING



// Preprocessor printing functions
#ifdef VERBOSE
#    define PRINT(outputStr) std::cout << "\n************************************************************\n" << "Making: " << (outputStr) << "\n" << "******\
******************************************************\n\n";
#    define SUBPRINT(outputStr) std::cout << "\t" << outputStr << "\n";
#    define QUIT(outputStr) std::cout << "\n\n\n" << outputStr << "\n\n\n"; exit(0);
#else
#    define PRINT(outputStr)
#    define SUBPRINT(outputStr)
#    define QUIT(outputStr)
#endif




#endif
