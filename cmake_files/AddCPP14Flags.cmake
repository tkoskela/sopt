include(CheckCXXCompilerFlag)

# On older cmake versions + newer compilers,
# the given version of CheckCXXCompilerFlags does not quite work.
if(CMAKE_VERSION VERSION_LESS 2.8.9)
  macro (CHECK_CXX_COMPILER_FLAG _FLAG _RESULT)
     set(SAFE_CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS}")
     set(CMAKE_REQUIRED_DEFINITIONS "${_FLAG}")
     CHECK_CXX_SOURCE_COMPILES("int main() { return 0;}" ${_RESULT}
       # Some compilers do not fail with a bad flag
       FAIL_REGEX "command line option .* is valid for .* but not for C\\\\+\\\\+" # GNU
       FAIL_REGEX "unrecognized .*option"                     # GNU
       FAIL_REGEX "unknown .*option"                          # Clang
       FAIL_REGEX "ignoring unknown option"                   # MSVC
       FAIL_REGEX "warning D9002"                             # MSVC, any lang
       FAIL_REGEX "option.*not supported"                     # Intel
       FAIL_REGEX "invalid argument .*option"                 # Intel
       FAIL_REGEX "ignoring option .*argument required"       # Intel
       FAIL_REGEX "[Uu]nknown option"                         # HP
       FAIL_REGEX "[Ww]arning: [Oo]ption"                     # SunPro
       FAIL_REGEX "command option .* is not recognized"       # XL
       FAIL_REGEX "not supported in this configuration; ignored"       # AIX
       FAIL_REGEX "File with unknown suffix passed to linker" # PGI
       FAIL_REGEX "WARNING: unknown flag:"                    # Open64
       )
     set (CMAKE_REQUIRED_DEFINITIONS "${SAFE_CMAKE_REQUIRED_DEFINITIONS}")
  endmacro ()
endif(CMAKE_VERSION VERSION_LESS 2.8.9)

check_cxx_compiler_flag(-std=c++14 has_std_cpp14)
check_cxx_compiler_flag(-std=c++0x has_std_cpp0x)
if(MINGW)
  check_cxx_compiler_flag(-std=gnu++14 has_std_gnupp14)
  check_cxx_compiler_flag(-std=gnu++0x has_std_gnupp0x)
endif(MINGW)
unset(CXX14_FLAGS)
if(has_std_gnupp14)
    set(CXX14_FLAGS "${CXX14_FLAGS} -std=gnu++14")
elseif(has_std_gnupp0x)
    set(CXX14_FLAGS "${CXX14_FLAGS} -std=gnu++0x")
elseif(has_std_cpp14)
    set(CXX14_FLAGS "${CXX14_FLAGS} -std=c++14")
elseif(has_std_cpp0x)
    set(CXX14_FLAGS "${CXX14_FLAGS} -std=c++0x")
endif(has_std_gnupp14)

if(MSVC)
  set(MSWINDOBE TRUE)
  add_definitions(/EHsc)
  # Wd4251 stops MSCrapWare from issuing meaningless warnings. Seems Microsoft engineers don't grok
  # dynamic libraries yet. Or templates. Or both acting alone or together. In any case, issuing
  # warning sure is easier on them than fixing  their OS.
  # Unfortunately, it does disable warnings that may be of interest. Possibly.
  set(CXX14_FLAGS "${CXX14_FLAGS} /D_VARIADIC_MAX=10 /wd4251")
endif(MSVC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX14_FLAGS}")
set(PROJECT_USES_CPP14 True CACHE INTERNAL "Uses c++14.")