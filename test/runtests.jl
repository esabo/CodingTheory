# If you've activated the environment of the package, then
#    ]test
# if you've installed the package, then
#    ]test CodingTheory
#
# testing starts a new julia process, uses the package and all test deps, then runs this file
#
# you can do things like:

using Tests

@test 5 == 2^2 + 1