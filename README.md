
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

The  type  of  signed  permutations  is  `SPerm{T}`  where  `T<:Integer`, a `struct`  with one  field, a  `Vector{T}` which  holds the  image of `1:n`. Using a `T` siferrent than `Int` may possibly save space or time. If `T` is not  specified we  take it  to be  `Int16` since  this is a good compromise between speed and compactness.

SPerms  have methods `copy, hash, ==, isless`  (total order) so they can be keys in hashes or elements of sets; two `SPerms` are equal if they move the same points to the same images. For instance,

```julia-repl
julia> SPerm([-2,-1,-3])==SPerm([-2,-1,-3,4])
true
```

SPerms are considered as scalars for broadcasting.


<a target='_blank' href='https://github.com/jmichel7/SignedPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L1-L45' class='documenter-source'>source</a><br>

<a id='SPerms.SPerm' href='#SPerms.SPerm'>#</a>
**`SPerms.SPerm`** &mdash; *Type*.



`struct SPerm`

An  `SPerm` represents a signed permutation of `1:n`, that is a permutation of  the  set  `-n,…,-1,1,…,n`  which  preserves  the  pairs `(-i,i)`. It is implemented  by a `struct SPerm` which has  one field `d`, a vector holding the images of `1:n`.


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L57-L64' class='documenter-source'>source</a><br>

<a id='SPerms.SPerm-Tuple{Vararg{Integer}}' href='#SPerms.SPerm-Tuple{Vararg{Integer}}'>#</a>
**`SPerms.SPerm`** &mdash; *Method*.



SPerm{T}(x::Integer...)where T<:Integer

returns   a   signed   cycle.   For  instance  `SPerm{Int8}(1,-2,-1,2)`  and `SPerm({Int8}[-2,1])`  define  the  same  signed  permutation. If not given `{T}` is taken to be `{Int16}`.


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L69-L75' class='documenter-source'>source</a><br>

<a id='SPerms.@sperm_str-Tuple{String}' href='#SPerms.@sperm_str-Tuple{String}'>#</a>
**`SPerms.@sperm_str`** &mdash; *Macro*.



@sperm"..."

makes a  `SPerm`  from  a  string  specifying  signed cycles linke the REPL printing of an `SPerm`; an example is `sperm"(1,-2)(5,-6,7)(-4,9)"`


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L92-L97' class='documenter-source'>source</a><br>

<a id='PermGroups.Perms.Perm-Tuple{SPerm}' href='#PermGroups.Perms.Perm-Tuple{SPerm}'>#</a>
**`PermGroups.Perms.Perm`** &mdash; *Method*.



`Perm(p::SPerm)` returns the underlying `Perm` of an `SPerm`


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L137' class='documenter-source'>source</a><br>

<a id='SPerms.signs' href='#SPerms.signs'>#</a>
**`SPerms.signs`** &mdash; *Function*.



`signs(p::SPerm)` returns the underlying signs of an `SPerm`


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L140' class='documenter-source'>source</a><br>

<a id='PermGroups.Perms.permute-Tuple{AbstractVector, SPerm}' href='#PermGroups.Perms.permute-Tuple{AbstractVector, SPerm}'>#</a>
**`PermGroups.Perms.permute`** &mdash; *Method*.



`permute(l::AbstractVector,p::SPerm)`

returns `l` permuted by `p`, a vector `r` such that `r[abs(i^p)]=l[i]*sign(i^p)`.

```julia-repl
julia> p=SPerm([-2,-1,-3])
SPerm{Int64}: (1,-2)(3,-3)

julia> permute([20,30,40],p)
3-element Vector{Int64}:
 -30
 -20
 -40
```

`permute`  can also act on  matrices with a keyword  `dims`. If `dims=1` it permutes  the lines, if `dims=2` the  columns and for `dims=(1,2)` both. If `P=Matrix(p)`  and  `iP=Matrix(inv(p))`  then  `permute(m,p;dims=1)==iP*m`, `permute(m,p;dims=2)==m*P`  and `permute(m,p;dims=(1,2))==iP*m*P`. Finally, the  form  `permute(m,p1,p2)`  permutes  the  lines  of `m` by `p1` and the columns by `p2`.


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L311-L333' class='documenter-source'>source</a><br>

<a id='PermGroups.Perms.orbit-Tuple{SPerm, Integer}' href='#PermGroups.Perms.orbit-Tuple{SPerm, Integer}'>#</a>
**`PermGroups.Perms.orbit`** &mdash; *Method*.



`orbit(a::SPerm,i::Integer)` returns the orbit of `a` on `i`.


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L184-L186' class='documenter-source'>source</a><br>

<a id='PermGroups.Perms.order-Tuple{SPerm}' href='#PermGroups.Perms.order-Tuple{SPerm}'>#</a>
**`PermGroups.Perms.order`** &mdash; *Method*.



`order(a::SPerm)` is the order of the signed permutation `a`.


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L253-L255' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L199-L210' class='documenter-source'>source</a><br>

<a id='PermGroups.Perms.cycletype-Tuple{SPerm}' href='#PermGroups.Perms.cycletype-Tuple{SPerm}'>#</a>
**`PermGroups.Perms.cycletype`** &mdash; *Method*.



`cycletype(p::SPerm,n=length(p.d))` pair  of  partitions  parameterizing  the  conjugacy  class  of  `p` in the hyperoctaedral group `Bₙ`

```julia-repl
julia> cycletype(SPerm(1,-1),2)
2-element Vector{Vector{Int64}}:
 [1]
 [1]
```


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L223-L233' class='documenter-source'>source</a><br>

<a id='Base.Matrix-Tuple{SPerm}' href='#Base.Matrix-Tuple{SPerm}'>#</a>
**`Base.Matrix`** &mdash; *Method*.



`Matrix(a::SPerm,n=length(a.d))`  permutation matrix for  `a` (operating on `n` points)

if `m=Matrix(a)` then `permutedims(v)*m==permute(v,a)`. Also `Diagonal(signs(a))*Matrix(Perm(a))==Matrix(a)`.

```julia-repl
julia> Matrix(SPerm([-2,-1,-3]))
3×3 Matrix{Int64}:
  0  -1   0
 -1   0   0
  0   0  -1
```


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L393-L406' class='documenter-source'>source</a><br>

<a id='SPerms.SPerm-Tuple{AbstractMatrix{<:Integer}}' href='#SPerms.SPerm-Tuple{AbstractMatrix{<:Integer}}'>#</a>
**`SPerms.SPerm`** &mdash; *Method*.



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


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L416-L431' class='documenter-source'>source</a><br>

<a id='SPerms.SPerm-Tuple{AbstractVector, AbstractVector}' href='#SPerms.SPerm-Tuple{AbstractVector, AbstractVector}'>#</a>
**`SPerms.SPerm`** &mdash; *Method*.



`SPerm{T}(l::AbstractVector,l1::AbstractVector)`

return  a `SPerm{T}` `p`  such that `permute(l,p)==l1`  if such `p` exists; returns  nothing otherwise.  If not  given `{T}`  is taken to be `{Int16}`. Needs the entries of `l` and `l1` to be sortable.

```julia-repl
julia> p=SPerm([20,30,40],[-40,-20,-30])
(1,-2,3,-1,2,-3)

julia> permute([20,30,40],p)
3-element Vector{Int64}:
 -40
 -20
 -30
```


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L363-L380' class='documenter-source'>source</a><br>

<a id='SPerms.sstab_onmats' href='#SPerms.sstab_onmats'>#</a>
**`SPerms.sstab_onmats`** &mdash; *Function*.



`sstab_onmats([G,]M[,l])`

If  `onmats(M,p)=permute(M,p,p)` (simultaneous  signed conjugation  of rows and  columns, or conjugating by the  matrix of the signed permutation `p`), and  the argument `G`  is given (which  should be an  `SPermGroup`) this is just  a fast implementation of `centralizer(G,M,onmats)`. If `G` is omitted it  is  taken  to  be  `hyperoctaedral_group(size(M,1))`.  The program uses sophisticated  algorithms, and can  handle matrices up  to 80×80. If `l` is given the return group should also centralize `l` (for the action ^)

```julia-repl
julia> n=[-1 -1 -1 -2 2 -2 -3 -3 -3; -1 -1 -1 -3 3 -3 -2 -2 -2; -1 -1 -1 -1 1 -1 -1 -1 -1; -2 -3 -1 -3 1 -2 -1 -3 -2; 2 3 1 1 -2 3 3 2 1; -2 -3 -1 -2 3 -1 -2 -1 -3; -3 -2 -1 -1 3 -2 -2 -3 -1; -3 -2 -1 -3 2 -1 -3 -1 -2; -3 -2 -1 -2 1 -3 -1 -2 -3];

julia> g=sstab_onmats(n)
Group([(1,6)(2,8)(5,-7),(1,8)(2,6)(4,9),(1,2)(4,9)(5,-7)(6,8),(1,-6)(2,-8)(3,-3)(4,-4)(5,7)(9,-9),(1,-8)(2,-6)(3,-3)(4,-9)(5,-5)(7,-7),(1,-2)(3,-3)(4,-9)(5,7)(6,-8)])

julia> length(g)
8
```


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L492-L512' class='documenter-source'>source</a><br>

<a id='SPerms.SPerm-Tuple{AbstractMatrix, AbstractMatrix}' href='#SPerms.SPerm-Tuple{AbstractMatrix, AbstractMatrix}'>#</a>
**`SPerms.SPerm`** &mdash; *Method*.



`SPerm(M::AbstractMatrix,N::AbstractMatrix;dims)`

returns a signed permutation `p` such that `permute(M,p;dims)==N` is such a `p`  exists,  and  `nothing`  otherwise.  If  `dims=(1,2)` then `M` and `N` should be symmetric matrices.

The  case `dims=(1,2)` routine is useful  to identify two objects which are isomorphic  but with different labelings. It  is used in Chevie to identify Lusztig  Fourier transform  matrices with  standard (classified)  data. The program  uses sophisticated algorithms, and can often handle matrices up to 80×80.

```julia-repl
julia> p=sperm"(1,-1)(2,5,3,-4,-2,-5,-3,4)(7,-9)";

julia> m=permute(n,p,p);

julia> p=SPerm(m,n;dims=(1,2))
(1,8,-1,-8)(2,-9,-7,-4,-3,5,-6)

julia> permute(m,p;dims=(1,2))==n
true
```


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L605-L629' class='documenter-source'>source</a><br>

<a id='SPerms.hyperoctaedral_group' href='#SPerms.hyperoctaedral_group'>#</a>
**`SPerms.hyperoctaedral_group`** &mdash; *Function*.



`hyperoctaedral_group(n)` the hyperoctaedral group on `1:n`


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L452' class='documenter-source'>source</a><br>

<a id='SPerms.SPerm_onmats' href='#SPerms.SPerm_onmats'>#</a>
**`SPerms.SPerm_onmats`** &mdash; *Function*.



`SPerm_onmats(M,N;extra=nothing)`

`M`  and `N` should be symmetric  matrices. `SPerm_onmats` returns a signed permutation `p` such that `onmats(M,p)=N` if such a permutation exists, and `nothing`  otherwise. If  in addition  vectors `extra=[m,n]` are given, the signed permutation `p` should also satisfy `permute(m,p)==n`.

Efficient version of `transporting_elt(hyperoctaedral_group(size(M,1)),M,N,onmats)`


<a target='_blank' href='https://github.com/jmichel7/SPerms.jl/blob/0f39c48960f0a5652792640e0263fd80a6d711e0/src/SPerms.jl#L534-L545' class='documenter-source'>source</a><br>

