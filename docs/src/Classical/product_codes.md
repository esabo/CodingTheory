# Product Codes

A `MatrixProductCode` combines constituent codes through a defining matrix: row
``i`` of the defining matrix says how codewords of the ``i``th constituent are
mixed into each output block. The ordinary product and tensor-product
constructions of two codes are in
[New codes from old](@ref new-codes-from-old-api).

```@autodocs
Modules = [CodingTheory]
Pages = ["Classical/MatrixProductCode.jl"]
Private = false
```
