# Linear Codes

## Constructors

Generic linear codes may be constructed in two ways: via a matrix or via a vector space object. If a vector space is used, the basis of the vector space is used as a generator matrix for the code. If the optional parameter `parity` is set to true, the input is considered a parity-check matrix instead of a generator matrix. At the moment, no convention is used for the zero code and an error is thrown for such imputs. Zero rows are automatically removed from the input but zero columns are not. See the [tutorials](link) for usage examples.
```@docs
LinearCode
```

## Attributes

Various getter/accessor functions are provided for accessing attributes about the codes. The user is strongly encouraged to use these functions and never to work with the underlying structs directly, as many functions rely on the information in the structs to be in a specific order and don't check if information has been updated.

```@docs
field
```



## Methods



