/** \file imiMacro.h
 *
 *  <small> <!--Copyright Information: -->
 *  \b Author: Jan Ehrhardt \n
 *  \b Copyright (C) 2010, Jan Ehrhardt, Institute of Medical Informatics,
 *     University of Luebeck, Germany. All rights reserved.\n
 *     Please cite the following paper:\n
 *     Ehrhardt, J., Säring, D., & Handels, H. (2007). Structure-preserving
 *     interpolation of temporal and spatial image sequences using an optical
 *     flow-based method. Methods of information in medicine, 46(03), 300-307.\n *
 *  </small>
 ****************************************************************************/
/*
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef __imiMacro_h
#define __imiMacro_h

// System includes:
#include <iostream>

namespace imi
{
/////////////////////////////////////////////////////////////
// Macros for printing messages.
/////////////////////////////////////////////////////////////
/** \brief Macro used in the imiLib for printing warnings. */
#define imiWARNING(x) std::cout<<"\nWarning: "<<x<<std::flush;

/** \brief Macro used in the imiLib for printing error messages. */
#define imiERROR(x) std::cout<<"\n\033[1;31mError: "<<x<<"\033[0m"<<std::flush;

/** \brief Macro used in the imiLib for printing debug messages. */
#ifdef IMI_NO_DEBUG
#define imiDEBUG(x)
#else
#define imiDEBUG(x) std::cerr<<"\nDebug: In " __FILE__ ", line " << __LINE__<<":"<<x<<std::flush;
#endif

/** \brief Macro used in the imiLib for printing info messages. */
#ifdef IMI_NO_VERBOSE
#define imiINFO(x)
#else
//#define imiINFO(x) std::cout << "\n" << (imi::TIME_INITIALIZED ? imi::GET_MONOTONIC_TIME() : "Info: ") << x << std::flush;
#define imiINFO(x) std::cout << "\nInfo: " << x << std::flush;
#endif

/** \brief Macro used in the imiLib for printing info messages, depending on the
 *  global debug level. */
#ifdef IMI_NO_DEBUG
#define imiDEBUGINFO(level,msg)
#else
#define imiDEBUGINFO(level,msg)  \
  { if( imiImageFunctions::CheckGlobalDebugLevel( level ) ) \
    { \
      std::cout << "\n" << "Info: " << msg << std::flush; \
    } \
  }
#endif

} // namespace imi

#endif
