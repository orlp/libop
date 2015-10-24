libop
=====

This is my (Orson Peters) personal C++ library. It contains a wide array of features that I tend to
use in projects. There is a large overlap between boost and my library, as well as future versions
of the standard library. For example I have a C++11 port of `std::optional` as `op::optional`. As
time progresses and standard library support improves I will remove these.

The library is released under the permissive zlib license.

The entire library is header-only, so installing is a simple
`git clone https://github.com/orlp/libop`. Updating the library is a simple `git pull`.

There is no real documentation to speak of, however the headers in `bits/` contain a short
description of each function / class.

