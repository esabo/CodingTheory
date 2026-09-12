# [Subsystem Codes and Shared Accessors](@id quantum-subsystem-api)

A subsystem code is a stabilizer code in which some logical qubits are
designated as gauge qubits and left unprotected. Because a stabilizer code is
the special case with no gauge qubits, this file supplies the accessors used by
every quantum code in the library, which is why the general accessors appear
here rather than on the [core API page](@ref quantum-code-api).

The distinction that runs through this page is between the stabilizer group,
the gauge group, and the logical operators. *Bare* logical operators commute
with the entire gauge group; *dressed* logical operators need only commute with
the stabilizer group, so they may be multiplied by gauge operators and can have
lower weight. Functions come in bare and dressed forms wherever the two differ,
and conflating them will give the wrong distance.

Operators are returned either as matrices in symplectic ``[X \mid Z]`` form or
as vectors of operator pairs, depending on the function; the docstrings say
which. Signs are tracked separately through the character vector.

`promote_gauges_to_logical` and `promote_logicals_to_gauge` move qubits between
the two roles, which is the usual way to trade protected qubits for a lower
measurement weight. The `!` forms modify the code in place.

```@autodocs
Modules = [CodingTheory]
Pages = ["Quantum/subsystem_code.jl"]
Private = false
```
