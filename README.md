
<a id='Signed-permutations'></a>

<a id='Signed-permutations-1'></a>

# Signed permutations

<a id='SignedPerms' href='#SignedPerms'>#</a>
**`SignedPerms`** &mdash; *Module*.



A  signed permutation of `1:n` is  a permutation of the set `-n,…,-1,1,…,n` which  preserves the  pairs `(-i,i)`.  It is  represented internally as the images of `1:n`. It is printed as a product of signed cycles.

**Examples**

```julia-repl
julia> SPerm([-2,-1,-3])
SPerm{Int64}: (1,-2)(3,-3)
```

The group of signed permutations of `1:n` is called the hyperoctaedral group.

```
julia> W=hyperoctaedral_group(2)
Group([(1,-1),(1,2)])

julia> elements(W)
8-element Vector{SPerm{Int8}}:
 ()
 (1,-1)
 (1,2)
 (1,-2,-1,2)
 (1,2,-1,-2)
 (1,-2)
 (2,-2)
 (1,-1)(2,-2)
```

A  motivation for my use of signed  permutations is to find if two matrices differ  only by a simultaneous signed permutation of lines and columns. See the example below with `SPerm(m,n;dims=(1,2))`.

The  type  of  signed  permutations  is  `SPerm{T}`  where  `T<:Integer`, a `struct`  with one  field, a  `Vector{T}` which  holds the  image of `1:n`. Using a `T` different than `Int` may possibly save space or time. If `T` is not  specified we  take it  to be  `Int16` since  this is a good compromise between speed and compactness.

`SPerm`s  have methods `copy, hash, ==, isless`  (total order) so they can be keys in hashes or elements of sets; two `SPerm`s are equal if they move the same points to the same images. For instance,

```julia-repl
julia> SPerm([-2,-1,-3])==SPerm([-2,-1,-3,4])
true
```

`SPerm`s are considered as scalars for broadcasting.


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L1-L45' class='documenter-source'>source</a><br>

<a id='SignedPerms.SPerm' href='#SignedPerms.SPerm'>#</a>
**`SignedPerms.SPerm`** &mdash; *Type*.



`struct SPerm`

An  `SPerm` represents a signed permutation of `1:n`, that is a permutation of  the  set  `-n,…,-1,1,…,n`  which  preserves  the  pairs `(-i,i)`. It is implemented  by a `struct SPerm` which has  one field `d`, a vector holding the images of `1:n`. It is printed as its cycle decomposition.


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L56-L63' class='documenter-source'>source</a><br>

<a id='SignedPerms.SPerm-Tuple{Vararg{Integer}}' href='#SignedPerms.SPerm-Tuple{Vararg{Integer}}'>#</a>
**`SignedPerms.SPerm`** &mdash; *Method*.



SPerm{T}(x::Integer...)where T<:Integer

returns as a `SPerm{T}` a signed cycle. For instance `SPerm{Int8}(1,-2,-1,2)`  and `SPerm({Int8}[-2,1])` define  the same signed permutation. If not given `{T}` is taken to be `{Int16}`.


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L70-L76' class='documenter-source'>source</a><br>

<a id='SignedPerms.@sperm_str-Tuple{String}' href='#SignedPerms.@sperm_str-Tuple{String}'>#</a>
**`SignedPerms.@sperm_str`** &mdash; *Macro*.



@sperm"..."

makes a  `SPerm`  from  a  string  specifying  signed cycles linke the REPL printing of an `SPerm`; an example is `sperm"(1,-2)(5,-6,7)(-4,9)"`


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L91-L96' class='documenter-source'>source</a><br>

<a id='PermGroups.Perms.Perm-Tuple{SPerm}' href='#PermGroups.Perms.Perm-Tuple{SPerm}'>#</a>
**`PermGroups.Perms.Perm`** &mdash; *Method*.



`Perm(p::SPerm)` returns the underlying `Perm` of an `SPerm`


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L136' class='documenter-source'>source</a><br>

<a id='SignedPerms.signs' href='#SignedPerms.signs'>#</a>
**`SignedPerms.signs`** &mdash; *Function*.



`signs(p::SPerm)` returns the underlying signs of an `SPerm`


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L139' class='documenter-source'>source</a><br>

<a id='PermGroups.Perms.invpermute-Tuple{AbstractVector, SPerm}' href='#PermGroups.Perms.invpermute-Tuple{AbstractVector, SPerm}'>#</a>
**`PermGroups.Perms.invpermute`** &mdash; *Method*.



`invpermute(l::AbstractVector,p::SPerm)`

returns `l` invpermuted by `p`, a vector `r` such that `r[abs(i^p)]=l[i]*sign(i^p)`.

```julia-repl
julia> p=SPerm([-2,-1,-3])
SPerm{Int64}: (1,-2)(3,-3)

julia> invpermute([20,30,40],p)
3-element Vector{Int64}:
 -30
 -20
 -40
```

`invpermute` can also act on matrices with a keyword `dims`. If `dims=1` it invpermutes  the lines, if `dims=2` the  columns and for `dims=(1,2)` both. If `P=Matrix(p)` and `iP=Matrix(inv(p))` then `invpermute(m,p;dims=1)==iP*m`,      `invpermute(m,p;dims=2)==m*P`      and `invpermute(m,p;dims=(1,2))==iP*m*P`. Finally, the form `invpermute(m,p1,p2)`  invpermutes the lines of `m` by `p1` and the columns by `p2`.


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L310-L334' class='documenter-source'>source</a><br>

<a id='PermGroups.Perms.orbit-Tuple{SPerm, Integer}' href='#PermGroups.Perms.orbit-Tuple{SPerm, Integer}'>#</a>
**`PermGroups.Perms.orbit`** &mdash; *Method*.



`orbit(p::SPerm,i::Integer)` returns the orbit of `p` on `i`.


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L189-L191' class='documenter-source'>source</a><br>

<a id='PermGroups.Perms.order-Tuple{SPerm}' href='#PermGroups.Perms.order-Tuple{SPerm}'>#</a>
**`PermGroups.Perms.order`** &mdash; *Method*.



`order(a::SPerm)` is the order of the signed permutation `a`.


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L251-L253' class='documenter-source'>source</a><br>

<a id='PermGroups.Perms.last_moved-Tuple{SPerm}' href='#PermGroups.Perms.last_moved-Tuple{SPerm}'>#</a>
**`PermGroups.Perms.last_moved`** &mdash; *Method*.



`last_moved(a::SPerm)` is the largest integer moved by a


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L183' class='documenter-source'>source</a><br>

<a id='PermGroups.Perms.cycles-Tuple{SPerm}' href='#PermGroups.Perms.cycles-Tuple{SPerm}'>#</a>
**`PermGroups.Perms.cycles`** &mdash; *Method*.



`cycles(p::SPerm)` the non-trivial cycles of `p`.

Two cycles which differ only by sign are returned once only.

```julia-repl
julia> cycles(SPerm(-1,2)*SPerm(3,-3)*SPerm(4,5,-4,-5))
3-element Vector{Vector{Int16}}:
 [1, -2]
 [3, -3]
 [4, 5, -4, -5]
```


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L202-L213' class='documenter-source'>source</a><br>

<a id='PermGroups.Perms.cycletype-Tuple{SPerm}' href='#PermGroups.Perms.cycletype-Tuple{SPerm}'>#</a>
**`PermGroups.Perms.cycletype`** &mdash; *Method*.



`cycletype(p::SPerm,n=last_moved(p))` pair  of  partitions  parameterizing  the  conjugacy  class  of  `p` in the hyperoctaedral group `Bₙ`

```julia-repl
julia> cycletype(SPerm(1,-1),2)
2-element Vector{Vector{Int64}}:
 [1]
 [1]
```


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L226-L236' class='documenter-source'>source</a><br>

<a id='Base.Matrix-Tuple{SPerm}' href='#Base.Matrix-Tuple{SPerm}'>#</a>
**`Base.Matrix`** &mdash; *Method*.



`Matrix(p::SPerm,n=last_moved(p))`  permutation matrix for `p` operating on `n` points.

For a vector `v`, we have `permutedims(v)*Matrix(p)==invpermute(v,p)`. Also `Diagonal(signs(p))*Matrix(Perm(p))==Matrix(p)`.

```julia-repl
julia> Matrix(SPerm([-2,-1,-3]))
3×3 Matrix{Int64}:
  0  -1   0
 -1   0   0
  0   0  -1
```


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L400-L413' class='documenter-source'>source</a><br>

<a id='SignedPerms.SPerm-Tuple{AbstractMatrix{<:Integer}}' href='#SignedPerms.SPerm-Tuple{AbstractMatrix{<:Integer}}'>#</a>
**`SignedPerms.SPerm`** &mdash; *Method*.



`SPerm{T}(m::AbstractMatrix)`  If  `m`  is  a  signed  permutation  matrix, returns  the corresponding signed permutation of  type `T`. If omitted, `T` is taken to be `Int16`.

```julia-repl
julia> m=[0 -1 0;-1 0 0;0 0 -1]
3×3 Matrix{Int64}:
  0  -1   0
 -1   0   0
  0   0  -1

julia> SPerm(m)
(1,-2)(3,-3)
```


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L422-L437' class='documenter-source'>source</a><br>

<a id='SignedPerms.SPerm-Tuple{AbstractVector, AbstractVector}' href='#SignedPerms.SPerm-Tuple{AbstractVector, AbstractVector}'>#</a>
**`SignedPerms.SPerm`** &mdash; *Method*.



`SPerm{T}(l::AbstractVector,l1::AbstractVector)`

return  a `SPerm{T}` `p`  such that `invpermute(l1,p)==l`  if such `p` exists; returns  nothing otherwise.  If not  given `{T}`  is taken to be `{Int16}`. Needs the entries of `l` and `l1` to be sortable.

```julia-repl
julia> p=SPerm([20,30,40],[-40,-20,-30])
(1,-3,2,-1,3,-2)

julia> invpermute([-40,-20,-30],p)
3-element Vector{Int64}:
 20
 30
 40
```


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L370-L387' class='documenter-source'>source</a><br>

<a id='PermGroups.Perms.onmats-Tuple{AbstractMatrix, SPerm}' href='#PermGroups.Perms.onmats-Tuple{AbstractMatrix, SPerm}'>#</a>
**`PermGroups.Perms.onmats`** &mdash; *Method*.



`onmats(m::AbstractMatrix,g::SPerm)` synonym for `invpermute(m,g;dims=(1,2))` or `invpermute(m,g,g)`.


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L344-L347' class='documenter-source'>source</a><br>

<a id='SignedPerms.sstab_onmats' href='#SignedPerms.sstab_onmats'>#</a>
**`SignedPerms.sstab_onmats`** &mdash; *Function*.



`sstab_onmats([G,]M[,l])`

if  the argument `G`  is given (which  should be an  `SPermGroup`) this is just  a fast implementation of `centralizer(G,M,onmats)`. If `G` is omitted it  is  taken  to  be  `hyperoctaedral_group(size(M,1))`.  The program uses sophisticated  algorithms, and can  handle matrices up  to 80×80. If `l` is given the return group should also centralize `l` (for the action ^)

```julia-repl
julia> n=[-1 -1 -1 -2 2 -2 -3 -3 -3; -1 -1 -1 -3 3 -3 -2 -2 -2; -1 -1 -1 -1 1 -1 -1 -1 -1; -2 -3 -1 -3 1 -2 -1 -3 -2; 2 3 1 1 -2 3 3 2 1; -2 -3 -1 -2 3 -1 -2 -1 -3; -3 -2 -1 -1 3 -2 -2 -3 -1; -3 -2 -1 -3 2 -1 -3 -1 -2; -3 -2 -1 -2 1 -3 -1 -2 -3];

julia> length(sstab_onmats(n))
8
```


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L497-L512' class='documenter-source'>source</a><br>

<a id='SignedPerms.SPerm-Tuple{AbstractMatrix, AbstractMatrix}' href='#SignedPerms.SPerm-Tuple{AbstractMatrix, AbstractMatrix}'>#</a>
**`SignedPerms.SPerm`** &mdash; *Method*.



`SPerm(M::AbstractMatrix,N::AbstractMatrix;dims)`

returns  a signed  permutation `p`  such that  `invpermute(N,p;dims)==M` if such  a `p` exists,  and `nothing` otherwise.  If `dims=(1,2)` then `M` and `N` should be symmetric matrices.

The  case `dims=(1,2)` routine is useful  to identify two objects which are isomorphic  but with different labelings. It  is used in Chevie to identify Lusztig  Fourier transform  matrices with  standard (classified)  data. The program  uses sophisticated algorithms, and can often handle matrices up to 80×80.

```julia-repl
julia> p=sperm"(1,-1)(2,5,3,-4,-2,-5,-3,4)(7,-9)";

julia> m=onmats(n,p);

julia> onmats(n,SPerm(m,n;dims=(1,2)))==m
true
```


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L605-L626' class='documenter-source'>source</a><br>

<a id='SignedPerms.hyperoctaedral_group' href='#SignedPerms.hyperoctaedral_group'>#</a>
**`SignedPerms.hyperoctaedral_group`** &mdash; *Function*.



`hyperoctaedral_group(n)` the hyperoctaedral group on `1:n`


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L457' class='documenter-source'>source</a><br>

<a id='SignedPerms.SPerm_onmats' href='#SignedPerms.SPerm_onmats'>#</a>
**`SignedPerms.SPerm_onmats`** &mdash; *Function*.



`SPerm_onmats(M,N;extra=nothing)`

`M`  and `N` should be symmetric  matrices. `SPerm_onmats` returns a signed permutation `p` such that `onmats(N,p)=M` if such a permutation exists, and `nothing`  otherwise. If  in addition  vectors `extra=[m,n]` are given, the signed permutation `p` should also satisfy `invpermute(n,p)==m`.

Efficient version of `transporting_elt(hyperoctaedral_group(size(M,1)),M,N,onmats)`


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/d19d81a88f43ee8c012ad77d52fa336b169ff949/src/SignedPerms.jl#L534-L545' class='documenter-source'>source</a><br>

