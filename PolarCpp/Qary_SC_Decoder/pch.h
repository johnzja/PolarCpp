// pch.h: In fact, pre-compilation is disabled.

#ifndef PCH_H
#define PCH_H

#ifndef ON_LINUX
    #include "framework.h"
#endif

#include "mex.h"
#include "../PolarCpp/GF.h"
#include "../PolarCpp/SC.h"
#include "../PolarCpp/SCList.h"

#endif //PCH_H
