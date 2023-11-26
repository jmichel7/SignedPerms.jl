# auto-generated tests from julia-repl docstrings
using Test, SignedPerms, PermGroups
function mytest(file::String,cmd::String,man::String)
  println(file," ",cmd)
  exec=repr(MIME("text/plain"),eval(Meta.parse(cmd)),context=:limit=>true)
  if endswith(cmd,";") exec="nothing" 
  else exec=replace(exec,r"\s*$"m=>"")
       exec=replace(exec,r"\s*$"s=>"")
  end
  if exec!=man 
    i=1
    while i<=lastindex(exec) && i<=lastindex(man) && exec[i]==man[i]
      i=nextind(exec,i)
    end
    print("exec=$(repr(exec[i:end]))\nmanl=$(repr(man[i:end]))\n")
  end
  exec==man
end
@testset "SignedPerms.jl" begin
@test mytest("SignedPerms.jl","SPerm([-2,-1,-3])","SPerm{Int64}: (1,-2)(3,-3)")
@test mytest("SignedPerms.jl","SPerm([-2,-1,-3])==SPerm([-2,-1,-3,4])","true")
@test mytest("SignedPerms.jl","cycles(SPerm(-1,2)*SPerm(3,-3)*SPerm(4,5,-4,-5))","3-element Vector{Vector{Int16}}:\n [1, -2]\n [3, -3]\n [4, 5, -4, -5]")
@test mytest("SignedPerms.jl","cycletype(SPerm(1,-1),2)","2-element Vector{Vector{Int64}}:\n [1]\n [1]")
@test mytest("SignedPerms.jl","p=SPerm([-2,-1,-3])","SPerm{Int64}: (1,-2)(3,-3)")
@test mytest("SignedPerms.jl","invpermute([20,30,40],p)","3-element Vector{Int64}:\n -30\n -20\n -40")
@test mytest("SignedPerms.jl","p=SPerm([20,30,40],[-40,-20,-30])","(1,-3,2,-1,3,-2)")
@test mytest("SignedPerms.jl","invpermute([-40,-20,-30],p)","3-element Vector{Int64}:\n 20\n 30\n 40")
@test mytest("SignedPerms.jl","Matrix(SPerm([-2,-1,-3]))","3×3 Matrix{Int64}:\n  0  -1   0\n -1   0   0\n  0   0  -1")
@test mytest("SignedPerms.jl","m=[0 -1 0;-1 0 0;0 0 -1]","3×3 Matrix{Int64}:\n  0  -1   0\n -1   0   0\n  0   0  -1")
@test mytest("SignedPerms.jl","SPerm(m)","(1,-2)(3,-3)")
@test mytest("SignedPerms.jl","n=[-1 -1 -1 -2 2 -2 -3 -3 -3; -1 -1 -1 -3 3 -3 -2 -2 -2; -1 -1 -1 -1 1 -1 -1 -1 -1; -2 -3 -1 -3 1 -2 -1 -3 -2; 2 3 1 1 -2 3 3 2 1; -2 -3 -1 -2 3 -1 -2 -1 -3; -3 -2 -1 -1 3 -2 -2 -3 -1; -3 -2 -1 -3 2 -1 -3 -1 -2; -3 -2 -1 -2 1 -3 -1 -2 -3];","nothing")
@test mytest("SignedPerms.jl","length(sstab_onmats(n))","8")
@test mytest("SignedPerms.jl","p=sperm\"(1,-1)(2,5,3,-4,-2,-5,-3,4)(7,-9)\";","nothing")
@test mytest("SignedPerms.jl","m=onmats(n,p);","nothing")
@test mytest("SignedPerms.jl","onmats(n,SPerm(m,n;dims=(1,2)))==m","true")
end
