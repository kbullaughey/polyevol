Getting Started
===============

Try this once inside the 'src' directory:

    make
    cd examples
    ./run.sh

Tests
-----

I use [Google Test][gtest] for performing xUnit-style unit testing. I include version 1.6.0 in the repository, and running `make` from `src` will also build gtest itself as well as the tests. Tests are located in `tests/`, and can be run as follows:

    test/runner

[gtest]: http://code.google.com/p/googletest/

Requirements
------------

The above requires that you have g++ and R installed. And you'll need both the `gsl` and `grid` R libraries. I've only tried building the code on linux and Mac OS X. I distribute the code with Google Test, which I use as my unit test framework (similar to xUnit). 

Kevin Bullaughey

Mon, July 2 2012
