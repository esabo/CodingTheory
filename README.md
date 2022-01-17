# CodingTheory
A basic coding theory library for Julia.

## TODO:
- Whole library is based on row operations. Switch over to the transpose since Julia is column major to get a major speedup.
- Use @views
- Remove manual optimizations which might hinder the compiler, such as not doing something if the result is zero.
    Basic principals to use here:
    1. Branches are expensive, adds are not
    2. If assignments are unconditional, the compiler and processor microcode can optimize them
    3. If it thinks the old value may still be needed then it will have to keep it around in a register
- Make Julia package
- Finish documentation
- Do unit testing
