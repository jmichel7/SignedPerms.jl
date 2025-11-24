"""
A  signed permutation of `1:n` is  a permutation of the set `-n,…,-1,1,…,n`
which  preserves the  pairs `(-i,i)`.  It is  represented internally as the
images of `1:n`. It is printed as a product of signed cycles.

# Examples
```julia-repl
julia> SPerm([-2,-1,-3])
SPerm{Int64}: (1,-2)(3,-3)
```
The group of signed permutations of `1:n` is called the hyperoctaedral group.
```julia-rep1
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
A  motivation for my use of signed  permutations is to find if two matrices
differ  only by a simultaneous signed permutation of lines and columns. See
the example below with `SPerm(m,n;dims=(1,2))`.

The type of signed permutations is [`SPerm`](@ref)`{T}` where `T<:Integer`,
a  `struct` with one field,  a `Vector{T}` which holds  the image of `1:n`.
Using a `T` different than `Int` may possibly save space or time. If `T` is
not  specified we  take it  to be  `$Idef` since  this is a good compromise
between speed and compactness.

`SPerm`s  have methods `copy, hash, ==, isless`  (total order) so they can be
keys in hashes or elements of sets; two `SPerm`s are equal if they move the
same points to the same images. For instance,
```julia-repl
julia> SPerm([-2,-1,-3])==SPerm([-2,-1,-3,4])
true
```
`SPerm`s are considered as scalars for broadcasting.

For more information, look at

[`SPerm`](@ref),
[`signs`](@ref),
[`invpermute`](@ref),
[`orbit`](@ref),
[`order`](@ref),
[`last_moved`](@ref),
[`cycles`](@ref),
[`cycletype`](@ref),
[`Matrix`](@ref),
[`onmats`](@ref),
[`sstab_onmats`](@ref),
[`hyperoctaedral_group`](@ref),
[`SPerm_onmats`](@ref).
"""
module SignedPerms
using PermGroups
using Combinat: tally, collectby
export SPerm, sstab_onmats, @sperm_str, signs, SPermGroup, hyperoctaedral_group

const info=Ref(true)
function InfoChevie(a...)
  if info[] print(a...) end
end

"""
`struct SPerm`

An  `SPerm` represents a signed permutation of `1:n`, that is a permutation
of  the  set  `-n,…,-1,1,…,n`  which  preserves  the  pairs `(-i,i)`. It is
implemented  by a `struct SPerm` which has  one field `d`, a vector holding
the images of `1:n`. It is printed as its cycle decomposition.
"""
struct SPerm{T<:Integer}
   d::Vector{T}
end

const Idef=Int16 # default type T for SPerm{T}

"""
SPerm{T}(x::Integer...)where T<:Integer

returns as a `SPerm{T}` a signed cycle. For instance
`SPerm{Int8}(1,-2,-1,2)`  and `SPerm({Int8}[-2,1])` define  the same signed
permutation. If not given `{T}` is taken to be `{$Idef}`.
"""
function SPerm{T}(x::Integer...)where T<:Integer
  if isempty(x) return SPerm(T[]) end
  d=T.(1:max(abs.(x)...))
  for i in 1:length(x)-1
   d[abs(x[i])]=sign(x[i])*x[i+1]
  end
  if length(x)==1 d[abs(x[1])]=x[1]
  else d[abs(x[end])]=sign(x[end])*x[1]
  end
  SPerm(d)
end

SPerm(x::Integer...)=SPerm{Idef}(x...)

"""
   @sperm"..."

makes a  `SPerm`  from  a  string  specifying  signed cycles linke the REPL
printing of an `SPerm`; an example is `sperm"(1,-2)(5,-6,7)(-4,9)"`
"""
macro sperm_str(s::String)
  start=1
  res=SPerm()
  while true
    m=match(r"\((\s*-?\d+\s*,)+\s*-?\d+\)",s[start:end])
    if isnothing(m) break end
    start+=m.match.ncodeunits
    res*=SPerm(Meta.parse(m.match).args...)
  end
  res::SPerm
end

Base.convert(::Type{SPerm{T}},p::SPerm{T1}) where {T,T1}=T==T1 ? p : SPerm(T.(p.d))

@GapObj struct SPermGroup{T}<:Group{SPerm{T}}
  gens::Vector{SPerm{T}}
  one::SPerm{T}
end

Base.one(G::SPermGroup)=G.one

function Base.show(io::IO,G::SPermGroup)
  print(io,"Group([")
  join(io,gens(G),',')
  print(io,"])")
end

function Groups.Group(a::AbstractVector{SPerm{T}}) where T
  SPermGroup(filter(!isone,a),one(prod(a)),Dict{Symbol,Any}())
end

SPermGroup()=SPermGroup(SPerm{Idef}[],SPerm{Idef}(),Dict{Symbol,Any}())

Base.one(p::SPerm{T}) where T=SPerm(collect(T(1):T(length(p.d))))
Base.one(::Type{SPerm{T}}) where T=SPerm(T[])
Base.copy(p::SPerm)=SPerm(copy(p.d))

Perms.perm(a::SPerm)=a.d

"`Perm(p::SPerm)` returns the underlying `Perm` of an `SPerm`"
Perms.Perm(p::SPerm)=Perm(abs.(p.d))
SPerm(p::Perm)=SPerm(perm(p))
"`signs(p::SPerm)` returns the underlying signs of an `SPerm`"
signs(p::SPerm)=sign.(p.d)
# if N=onmats(M,p) then M==onmats(N^Diagonal(signs(p)),Perm(p))

# SPerms are scalars for broadcasting"
Base.broadcastable(p::SPerm)=Ref(p)

# hash is needed for using SPerms in Sets/Dicts
function Base.hash(a::SPerm, h::UInt)
  for (i,v) in pairs(a.d)
    if v!=i h=hash(v,h) end
  end
  h
end

function Base.promote_rule(a::Type{SPerm{T1}},b::Type{SPerm{T2}})where {T1,T2}
  SPerm{promote_type(T1,T2)}
end

extend!(a::SPerm,n::Integer)=if length(a.d)<n append!(a.d,length(a.d)+1:n) end

"""
 `promote_degree(a::SPerm, b::SPerm)` promotes `a` and `b` to the same type,
 then extends `a` and `b` to the same degree
"""
function promote_degree(a::SPerm,b::SPerm)
  a,b=promote(a,b)
  extend!(a,length(b.d))
  extend!(b,length(a.d))
  (a,b)
end

# total order is needed to use SPerms in sorted lists
function Base.isless(a::SPerm, b::SPerm)
  a,b=promote_degree(a,b)
  for (ai,bi) in zip(a.d,b.d) ai!=bi && return ai<bi end
  false
end

function Base.:(==)(a::SPerm, b::SPerm)
  a,b=promote_degree(a,b)
  a.d==b.d
end

" `last_moved(a::SPerm)` is the largest integer moved by a"
function Perms.last_moved(a::SPerm{T})where T
  @inbounds p=findlast(x->a.d[x]!=x,eachindex(a.d))
  isnothing(p) ? T(0) : T(p)
end

"""
`orbit(p::SPerm,i::Integer)` returns the orbit of `p` on `i`.
"""
function Perms.orbit(p::SPerm,i::Integer)
  res=[i]
  j=i
  while true
    j^=p
    if j==i return res end
    push!(res,j)
  end
end

"""
`cycles(p::SPerm)` the non-trivial cycles of `p`.

Two cycles which differ only by sign are returned once only.
```julia-repl
julia> cycles(SPerm(-1,2)*SPerm(3,-3)*SPerm(4,5,-4,-5))
3-element Vector{Vector{$Idef}}:
 [1, -2]
 [3, -3]
 [4, 5, -4, -5]
```
"""
function Perms.cycles(p::SPerm{T})where T
  orbs=Vector{T}[]
  to_visit=trues(length(p.d))
  for i in eachindex(to_visit)
    if !to_visit[i] continue end
    cyc=orbit(p,T(i))
@inbounds to_visit[abs.(cyc)].=false
    if length(cyc)>1 push!(orbs,cyc) end
  end
  orbs
end

"""
`cycletype(p::SPerm,n=last_moved(p))`
pair  of  partitions  parameterizing  the  conjugacy  class  of  `p` in the
hyperoctaedral group `Bₙ`
```julia-repl
julia> cycletype(SPerm(1,-1),2)
2-element Vector{Vector{Int64}}:
 [1]
 [1]
```
"""
function Perms.cycletype(p::SPerm,n=last_moved(p))
  res=[Int[],Int[]]
  to_visit=trues(n)
  for i in 1:n
    if !to_visit[i] continue end
    cyc=orbit(p, i)
    if -i in cyc push!(res[2], div(length(cyc),2))
    else push!(res[1], length(cyc))
    end
@inbounds to_visit[abs.(cyc)].=false
  end
  sort!.(res,rev=true)
end

"""
`order(a::SPerm)` is the order of the signed permutation `a`.
"""
PermGroups.order(a::SPerm) = lcm(length.(cycles(a)))

function Base.show(io::IO, a::SPerm)
  hasdecor=get(io,:TeX,false)||get(io,:limit,false)
  if !hasdecor print(io,"sperm\"") end
  cyc=cycles(a)
  if isempty(cyc) print(io,"()")
  else for c in cyc print(io,"(",join(c,","),")") end
  end
  if !hasdecor print(io,"\"") end
end

function Base.show(io::IO, ::MIME"text/plain", p::SPerm{T})where T
  if T!=Idef && !haskey(io,:typeinfo) print(io,typeof(p),": ") end
  show(io,p)
end

function Base.:*(a::SPerm, b::SPerm)
  a,b=promote_degree(a,b)
  r=similar(a.d)
  for (i,v) in pairs(a.d)
@inbounds if v<0 r[i]=-b.d[-v] else r[i]=b.d[v] end
  end
  SPerm(r)
end

function Base.inv(a::SPerm)
  r=similar(a.d)
  for (i,v) in pairs(a.d)
@inbounds if v<0 r[-v]=-i  else r[v]=i end
  end
  SPerm(r)
end

# less allocations than inv(a)*b
function Base.:\(a::SPerm, b::SPerm)
  a,b=promote_degree(a,b)
  r=similar(a.d)
  for (i,v) in pairs(a.d) 
@inbounds if v<0 r[-v]=-b.d[i] else r[v]=b.d[i] end
  end
  SPerm(r)
end

Base.:/(a::SPerm,b::SPerm)=a*inv(b)
Base.:^(a::SPerm, b::SPerm)=inv(b)*a*b
Base.:^(a::SPerm, b::Perm)=a^SPerm(b)

@inline function Base.:^(n::T, a::SPerm) where T<:Integer
  if abs(n)>length(a.d) return n end
@inbounds n<0 ? T(-a.d[-n]) : T(a.d[n])
end

Base.:^(a::SPerm, n::Integer)=n>=0 ? Base.power_by_squaring(a,n) :
                                     Base.power_by_squaring(inv(a),-n)

"""
`invpermute(l::AbstractVector,p::SPerm)`

returns `l` invpermuted by `p`, a vector `r` such that
`r[abs(i^p)]=l[i]*sign(i^p)`.

```julia-repl
julia> p=SPerm([-2,-1,-3])
SPerm{Int64}: (1,-2)(3,-3)

julia> invpermute([20,30,40],p)
3-element Vector{Int64}:
 -30
 -20
 -40
```

`invpermute` can also act on matrices with a keyword `dims`. If `dims=1` it
invpermutes  the lines, if `dims=2` the  columns and for `dims=(1,2)` both.
If `P=Matrix(p)` and `iP=Matrix(inv(p))` then
`invpermute(m,p;dims=1)==iP*m`,      `invpermute(m,p;dims=2)==m*P`      and
`invpermute(m,p;dims=(1,2))==iP*m*P`. Finally, the form
`invpermute(m,p1,p2)`  invpermutes the lines of `m` by `p1` and the columns
by `p2`.
"""
function Perms.invpermute(l::AbstractVector,a::SPerm)
  res=collect(l)
  for i in eachindex(l)
    v=i^a
    if v>0 res[v]=l[i] else res[-v]=-l[i] end
  end
  res
end

"""
`onmats(m::AbstractMatrix,g::SPerm)` synonym for `invpermute(m,g;dims=(1,2))`
or `invpermute(m,g,g)`.
"""
Perms.onmats(m::AbstractMatrix,g::SPerm)=invpermute(m,g,g)

function Perms.invpermute(m::AbstractMatrix,a::SPerm,b::SPerm)
  res=collect(m)
  for i in CartesianIndices(m)
    v1=getindex(i,1)^a
    v2=getindex(i,2)^b
    res[abs(v1),abs(v2)]=m[i]*sign(v1)*sign(v2)
  end
  res
end

function Perms.invpermute(m::AbstractMatrix,a::SPerm;dims=1)
  if dims==1 invpermute(m,a,SPerm())
  elseif dims==2 invpermute(m,SPerm(),a)
  elseif dims==(1,2) invpermute(m,a,a)
  end
end

# to find orbits under SPerms take abs
myabs(x)=x<-x ? -x : x

"""
`SPerm{T}(l::AbstractVector,l1::AbstractVector)`

return  a `SPerm{T}` `p`  such that `invpermute(l1,p)==l`  if such `p` exists;
returns  nothing otherwise.  If not  given `{T}`  is taken to be `{$Idef}`.
Needs the entries of `l` and `l1` to be sortable and have operation `-`.

```julia-repl
julia> p=SPerm([20,30,40],[-40,-20,-30])
(1,-3,2,-1,3,-2)

julia> invpermute([-40,-20,-30],p)
3-element Vector{Int64}:
 20
 30
 40
```
"""
function SPerm{T}(a::AbstractVector,b::AbstractVector)where T<:Integer
  p=Perm(myabs.(a),myabs.(b))
  if isnothing(p) return p end
  res=perm(p)
  for i in eachindex(a)
@inbounds if b[i]!=a[i^p] res[i]=-res[i] end
  end
  SPerm{T}(res)
end

SPerm(l::AbstractVector,l1::AbstractVector)=SPerm{Idef}(l,l1)

"""
`Matrix(p::SPerm,n=last_moved(p))`  permutation matrix for `p` operating on
`n` points.

For a vector `v`, we have `permutedims(v)*Matrix(p)==invpermute(v,p)`.
Also `Diagonal(signs(p))*Matrix(Perm(p))==Matrix(p)`.
```julia-repl
julia> Matrix(SPerm([-2,-1,-3]))
3×3 Matrix{Int64}:
  0  -1   0
 -1   0   0
  0   0  -1
```
"""
function Base.Matrix(a::SPerm,n=last_moved(a))
  res=zeros(Int,n,n)
  for i in 1:n res[i,abs(i^a)]=sign(i^a) end
  res
end

SPerm(m::AbstractMatrix)=SPerm{Idef}(m)

"""
`SPerm{T}(m::AbstractMatrix)`  If  `m`  is  a  signed  permutation  matrix,
returns  the corresponding signed permutation of  type `T`. If omitted, `T`
is taken to be `$Idef`.

```julia-repl
julia> m=[0 -1 0;-1 0 0;0 0 -1]
3×3 Matrix{Int64}:
  0  -1   0
 -1   0   0
  0   0  -1

julia> SPerm(m)
(1,-2)(3,-3)
```
"""
function SPerm{T}(m::AbstractMatrix) where T<:Integer
  n=size(m,1)
  if n!=size(m,2) error("matrix should be square") end
  res=fill(T(0),n)
  for i in 1:n
    if count(!iszero,@view m[i,:])!=1 error("not a monomial matrix") end
    nz=findfirst(!iszero,@view m[i,:])
    if m[i,nz]==1 res[i]=nz
    elseif m[i,nz]==-1 res[i]=-nz
    else error("not a signed permutation matrix")
    end
  end
  SPerm{T}(res)
end

# We have the property p=SPerm(Perm(p).d.*signs(p))

randSPerm(n)=SPerm(Idef.(rand((-1,1),n).*sortperm(rand(1:n,n))))

"""
`hyperoctaedral_group(n)` 

the hyperoctaedral group on `1:n`
```julia-repl
julia> W=hyperoctaedral_group(2)
Group([(1,-1),(1,2)])

julia> W(1,2)
SPerm{Int8}: (1,-2,-1,2)
```
"""
hyperoctaedral_group(n::Int)=
  Group(pushfirst!(map(i->SPerm{Int8}(i-1,i),2:n),SPerm{Int8}(1,-1)))

#--------------------- action on matrices -----------------------------------

# duplicate lines and cols of M so group(dup(...)) operates
function dup(M::AbstractMatrix)
  res=zeros(eltype(M),size(M).*2)
  for i in axes(M,1), j in axes(M,2)
    res[[2*i-1,2*i],[2*j-1,2*j]]=[M[i,j] -M[i,j];-M[i,j] M[i,j]]
  end
  res
end

dedup(M::AbstractMatrix)=M[1:2:size(M,1),1:2:size(M,2)]

# transform SPerm on -n:n to hyperoctaedral Perm acting on  1:2n
function dup(p::SPerm)
  res=empty(p.d)
  for i in p.d
    if i>0 push!(res,2i-1);push!(res,2i)
    else   push!(res,-2i);push!(res,-2i-1)
    end
  end
  Perm(res)
end

# transform hyperoctaedral Perm acting on  1:2n to SPerm on -n:n 
dedup(p::Perm{T}) where T=SPerm{T}(map(i->iseven(i) ? -div(i,2) : div(i+1,2),p.d[1:2:end-1]))

dup(g::SPermGroup)=Group(dup.(gens(g)))

Base.length(g::SPermGroup)=length(Group(dup.(gens(g))))

function invblocks(m,extra=nothing)
  if isnothing(extra) extra=zeros(Int,size(m,1)) end
  blk1=[collect(axes(m,1))]
  while true
    blk=blk1
    blk1=vcat(map(I->collectby(map(i->
            (tally(myabs.(m[i,I])),m[i,i],extra[i]),I),I), blk)...)
    if blk==blk1 return blk end
  end
end

"""
`sstab_onmats([G,]M[,l])`

if  the argument `G`  is given (which  should be an  `SPermGroup`) this is
just  a fast implementation of `centralizer(G,M,onmats)`. If `G` is omitted
it  is  taken  to  be  `hyperoctaedral_group(size(M,1))`.  The program uses
sophisticated  algorithms, and can  handle matrices up  to 80×80. If `l` is
given the return group should also centralize `l` (for the action ^)

```julia-repl
julia> n=[-1 -1 -1 -2 2 -2 -3 -3 -3; -1 -1 -1 -3 3 -3 -2 -2 -2; -1 -1 -1 -1 1 -1 -1 -1 -1; -2 -3 -1 -3 1 -2 -1 -3 -2; 2 3 1 1 -2 3 3 2 1; -2 -3 -1 -2 3 -1 -2 -1 -3; -3 -2 -1 -1 3 -2 -2 -3 -1; -3 -2 -1 -3 2 -1 -3 -1 -2; -3 -2 -1 -2 1 -3 -1 -2 -3]
9×9 Matrix{Int64}:
 -1  -1  -1  -2   2  -2  -3  -3  -3
 -1  -1  -1  -3   3  -3  -2  -2  -2
 -1  -1  -1  -1   1  -1  -1  -1  -1
 -2  -3  -1  -3   1  -2  -1  -3  -2
  2   3   1   1  -2   3   3   2   1
 -2  -3  -1  -2   3  -1  -2  -1  -3
 -3  -2  -1  -1   3  -2  -2  -3  -1
 -3  -2  -1  -3   2  -1  -3  -1  -2
 -3  -2  -1  -2   1  -3  -1  -2  -3

julia> G=sstab_onmats(n)
Group([(1,8)(2,6)(4,9),(1,6)(2,8)(5,-7),(1,-2)(3,-3)(4,-9)(5,7)(6,-8),(1,-6)(2,-8)(3,-3)(4,-4)(5,7)(9,-9),(1,2)(4,9)(5,-7)(6,8),(1,-8)(2,-6)(3,-3)(4,-9)(5,-5)(7,-7)])

julia> length(G)
8
```
"""
function sstab_onmats(M,extra=nothing)
  k=size(M,1)
  if M!=permutedims(M) error("M should be symmetric") end
  if isnothing(extra) extra=fill(1,size(M,1)) end
  blocks=sort(invblocks(M),by=length)
  gen=SPerm{Idef}[]
  I=Int[]
  for r in blocks
    if length(r)>5 InfoChevie("#IS Large Block:",r,"\n") end
    gr=stab_onmats(dup(hyperoctaedral_group(length(r))),dup(M[r,r]))
    p=SPerm(mappingPerm(1:length(r),r).d)
    append!(gen,map(x->dedup(x)^p,gens(gr)))
    append!(I,r)
    p=SPerm(mappingPerm(I,eachindex(I)).d)
    gen=gen.^p
    gen=dedup.(gens(stab_onmats(Group(dup.(gen)),dup(M[I,I]))))
    gen=gen.^inv(p)
  end
  return Group(gen)
end

"""
`SPerm_onmats(M,N;extra=nothing)`

`M`  and `N` should be symmetric  matrices. `SPerm_onmats` returns a signed
permutation `p` such that `onmats(N,p)=M` if such a permutation exists, and
`nothing`  otherwise. If  in addition  vectors `extra=[m,n]` are given, the
signed permutation `p` should also satisfy `invpermute(n,p)==m`.

Efficient version of
`transporting_elt(hyperoctaedral_group(size(M,1)),M,N,onmats)`

"""
function SPerm_onmats(M,N;extra=nothing)
  if M!=permutedims(M) error("M should be symmetric") end
  if N!=permutedims(N) error("N should be symmetric") end
  if isnothing(extra) extra=[fill(1,size(M,1)),fill(1,size(M,1))] end
  function ind(I,J)
    local iM,iN,p,n
    invM=map(i->(tally(myabs.(M[i,I])),M[i,i],extra[1][i]),I)
    invN=map(i->(tally(myabs.(N[i,J])),N[i,i],extra[2][i]),J)
    if tally(invM)!=tally(invN) InfoChevie("content differs");return end
    iM=collectby(invM,I)
    iN=collectby(invN,J)
    if length(iM)==1
      if length(I)>6 InfoChevie("large block:",length(I),"\n")
        p=transporting_elt(hyperoctaedral_group(length(I)),M[I,I],N[J,J],onmats,
            dist=(M,N)->count(i->M[i]!=N[i],eachindex(M)))
      else p=transporting_elt(hyperoctaedral_group(length(I)),M[I,I],N[J,J],onmats)
      end
      if isnothing(p) InfoChevie("could not match block\n");return end
      return [(I,J,p)]
    else p=map(ind,iM,iN)
      if nothing in p return  else return vcat(p...) end
    end
  end
  l=ind(axes(M,1),axes(N,1))
  if isnothing(l) return end
  I=Int[];J=Int[];g=SPermGroup();tr=SPerm()
  sort!(l)
  for r in l
#   @show r
    n=length(r[1])
#   q=mappingPerm(eachindex(I),eachindex(I))
    q=SPerm()
    p=SPerm(mappingPerm(1:n,(1:n).+length(I)))
    append!(I,r[1]);append!(J,r[2]);
#   @show I,J
    if !isone(comm(r[3]^p,tr^q)) error("noncomm") end
    tr=tr^q*r[3]^p
    h=gens(sstab_onmats(M[r[1],r[1]])).^p
#   @show h
    g=Group(vcat(gens(g).^q,h))
#   @show g
    e=transporting_elt(g,M[I,I],onmats(N[J,J],inv(tr)),onmats)
    if isnothing(e) 
#     @show g,I,J,inv(tr)
      return
    else
      if !isone(e^-1*e^tr) print("*** tr does not commute to e\n") end
      tr=e*tr
    end
#   @show I,J,e,tr
    g=stab_onmats(Group(dup.(gens(g))),dup(M[I,I]))
    g=Group(dedup.(gens(g)))
#   println(" #stab=",g)
  end
  # transporter of a ps from [1..Length(I)] to I
  trans=I->SPerm(mappingPerm(eachindex(I),I).d)
  trans(J)^-1*inv(tr)*trans(I)
end

"""
`SPerm(M::AbstractMatrix,N::AbstractMatrix;dims)`

returns  a signed  permutation `p`  such that  `invpermute(N,p;dims)==M` if
such  a `p` exists,  and `nothing` otherwise.  If `dims=(1,2)` then `M` and
`N` should be symmetric matrices.

The  case `dims=(1,2)` routine is useful  to identify two objects which are
isomorphic  but with different labelings. It  is used in Chevie to identify
Lusztig  Fourier transform  matrices with  standard (classified)  data. The
program  uses sophisticated algorithms, and can often handle matrices up to
80×80.

```julia-repl
julia> p=sperm"(1,-1)(2,5,3,-4,-2,-5,-3,4)(7,-9)";

julia> m=onmats(n,p);

julia> onmats(n,SPerm(m,n;dims=(1,2)))==m
true
```
"""
function SPerm{T}(m::AbstractMatrix,m1::AbstractMatrix;dims=1)where T<:Integer
  if     dims==1 SPerm{T}(collect(eachrow(m)),collect(eachrow(m1)))
  elseif dims==2 SPerm{T}(collect(eachcol(m)),collect(eachcol(m1)))
  elseif dims==(1,2) SPerm_onmats(m,m1)
  end
end

SPerm(m::AbstractMatrix,m1::AbstractMatrix;dims=1)=SPerm{Idef}(m,m1,dims=dims)

end
