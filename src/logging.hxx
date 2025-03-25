#ifndef __LOGGING_HXX__
#define __LOGGING_HXX__

#ifdef __VERBOSE_LOGGING__ // verbose logging, outputs a ton of information which may not be linear time to compute but useful for debugging on small graphs
#define __LOGGING__
#define V_LOG(a) std::cout << a; 
#else
#define V_LOG(a)
#endif

#ifdef __LOGGING__  // standard logging, outputs some additional information which is linear time to compute
#define __LIGHT_LOGGING__
#define N_LOG(a) std::cout << a;
#else
#define N_LOG(a)
#endif

#ifdef __LIGHT_LOGGING__ // light logging, outputs the bare minimum so you can verify things are working
#include <iostream>
#define L_LOG(a) std::cout << a;
#else
#define L_LOG(a)
#endif

#ifdef __DEBUG_LOGGING__ // debug logging, doesn't actually appear anywhere in the code but I sometimes slap some of these in there to debug and then it's easy to ctrl-f + delete all them after
						 // also flushes output in case of segfault so you can actually see it
#include <iostream>
#define D_LOG(a) std::cout << a; std::cout.flush();
#else
#define D_LOG(a)
#endif

#endif