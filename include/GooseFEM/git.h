/**
\file git.h
\brief Git information at install time.
*/

/**
Git commit hash.

Note that the actual value is filled by
-   CMake: at install time.
    ``CMakeLists.txt`` overwrite this file.
-   Python: at build time.
    ``setup.py`` defines this variable, this files is then never included on config.h.

If you install from conda the variables are filled at the server.
*/
#define GOOSEFEM_GIT_HASH "None"

/**
Git commit branch.

Note that the actual value is filled by
-   CMake: at install time.
    ``CMakeLists.txt`` overwrite this file.
-   Python: at build time.
    ``setup.py`` defines this variable, this files is then never included on config.h.

If you install from conda the variables are filled at the server.
*/
#define GOOSEFEM_GIT_BRANCH "None"
