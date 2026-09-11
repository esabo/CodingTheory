# Decoding LDPC Codes

## Message Passing
```@docs
SoftDecisionWorkspace
init_soft_workspace
load_soft_channel!
decode!
boxplus_exact
boxplus_minsum
boxplus_minsum_correction
csr_of
```

The soft-decision decoder selects sum-product and min-sum variants through
the `algorithm` keyword to `decode!`. Construct one reusable workspace per
thread; see [Message-passing Decoding](@ref message-passing-tutorial) for a
complete workflow.

## Linear Programming

```@docs
LP_decoder_LDPC
```

## Post-processing decoders

```@docs
OSDWorkspace
init_osd_workspace
osd_decode!
init_grand_workspace
grand_decode!
init_wbf_workspace
wbf_decode!
```