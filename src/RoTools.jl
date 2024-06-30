module RoTools

using CSV
using DataFrames
using Interpolations
using FilePathsBase

export BladeGeom, itr2d, hoverFM

include("geom_process.jl")
include("dust_pre_process.jl")
include("dust_post_process.jl")
include("flowuns_pre_process.jl")
include("flowuns_post_process.jl")

# 桨叶弦长和扭转分布规范化
function itr2d(x,y; extrapolation_method = Line()) # make convience LinearInterpolation of 2d array 
    res = Interpolations.LinearInterpolation(x,y,extrapolation_bc=extrapolation_method)
    return res 
end 

function itr2d(mat::Matrix; extrapolation_method = Line())
    x = mat[:,1]
    y = mat[:,2]
    res = Interpolations.LinearInterpolation(x,y,extrapolation_bc=extrapolation_method)
    return res 
end 

## 按照nrbe分段，重新绘制桨叶的几何外形（弦长分布和扭转分布）
## 输入：：初始桨叶几何形状的矩阵
## 输出：：按照nrbe分段数量重新内插而形成的新桨叶几何形状矩阵
function BladeGeom(mat::Matrix; nrbe=20)
    r = mat[:,1] # radius
    c = mat[:,2] # chord
    t = mat[:,3] # twist 

    itr_c = itr2d(r, c)
    itr_t = itr2d(r, t)

    rbe = range(r[1],r[end],nrbe)
    chord = itr_c.(rbe)
    twist = itr_t.(rbe)

    geom = hcat(rbe, chord, twist)

    return geom 
end     

# 用BSPL样条曲线来拟合一个标准的2列矩阵
# 函数`geom2bspl`：输入1个矩阵或者两个一维数组，返回`BSplines`插值类型的曲线/参数
any2flt(vec::Vector) = convert(Vector{Float64},vec) # 将数组转化为Float64的类型

function geom2bspl(mat::Matrix;norder=4,nctrp=3)
    basis = BSplines.BSplineBasis(norder,0:nctrp-1)
    spl = BSplines.interpolate(basis, any2flt(mat[:,1]), any2flt(mat[:,2]))
    bsplcoeffs = BSplines.coeffs(spl)
    return bsplcoeffs
end 

function geom2bspl(vec1::Vector,vec2::Vector;norder=4,nctrp=3)
    basis = BSplines.BSplineBasis(norder,0:nctrp-1)
    spl = BSplines.interpolate(basis, any2flt(vec1), any2flt(vec2))
    bsplcoeffs = BSplines.coeffs(spl)
    return bsplcoeffs
end 

# convert BSplines coffs to dust geom matrix 
function heli_blade_geom(opt_input,dustgeom)

    rbe = dustgeom.rbe
    bspl_chord = BSplines.BSplineBasis(4,range(rbe[1],rbe[end],3))
    bspl_twist = BSplines.BSplineBasis(3,range(rbe[1],rbe[end],2))
    cspl_chord = opt_input[1:5]
    cspl_twist = opt_input[6:8]

    cspl = BSplines.Spline(bspl_chord, cspl_chord)
    tspl = BSplines.Spline(bspl_twist, cspl_twist)

    chord = cspl.(rbe)
    twist = tspl.(rbe) 

    sweep = dustgeom.sweepback  
    dihed = dustgeom.dihedral
    n_seg = dustgeom.nseg 
    type_seg = dustgeom.typeseg 
    airfoil = dustgeom.aftable 
    geom = hcat(rbe,chord,twist,sweep,dihed,n_seg,type_seg,airfoil)
    dustgeom.rbe = rbe 
    dustgeom.chord = chord
    dustgeom.twist = twist 

    return geom, dustgeom 
end

# convert BSplines coffs to geom matrix 
function blade_geom_spl2std(spl_inp,rbe)
    bspl_chord = BSplines.BSplineBasis(4,range(rbe[1],rbe[end],3))
    bspl_twist = BSplines.BSplineBasis(3,range(rbe[1],rbe[end],3))
    cspl_chord = spl_inp[1:5]
    cspl_twist = spl_inp[6:9]

    cspl = BSplines.Spline(bspl_chord, cspl_chord)
    tspl = BSplines.Spline(bspl_twist, cspl_twist)

    chord = cspl.(rbe)
    twist = tspl.(rbe) 

    geom = hcat(collect(rbe),map(x->round(x,digits=4),chord), map(x->round(x,digits=4),twist))

    return geom 
end 

# convert geom matrix to BSplines coffs 
function blade_geom_std2spl(geom0::Matrix)
    geom = convert(Matrix{Float64},geom0[:,1:3])
    rbe = geom[:,1]
    chord = geom[:,2]
    twist = geom[:,3]
    bspl_chord = BSplines.BSplineBasis(4,range(rbe[1],rbe[end],3))
    bspl_twist = BSplines.BSplineBasis(3,range(rbe[1],rbe[end],3))
    spl_chord = BSplines.interpolate(bspl_chord,rbe,chord)
    spl_twist = BSplines.interpolate(bspl_twist,rbe,twist)

    coeffs_chord = BSplines.coeffs(spl_chord)
    coeffs_twist = BSplines.coeffs(spl_twist)

    return vcat(coeffs_chord, coeffs_twist)
end 

# 基于B样条曲线的桨叶几何外形DUST描述
function bspl_dust_geom(opt_input;afname="af01",nseg=1)

    rbe = 0.3:0.05:1.0
    nrow = length(rbe)
    bspl_chord = BSplines.BSplineBasis(4,0:2)
    bspl_twsit = BSplines.BSplineBasis(3,0:1)
    cspl_chord = opt_input[1:5]
    cspl_twist = opt_input[6:8]

    cspl = BSplines.Spline(bspl_chord, cspl_chord)
    tspl = BSplines.Spline(bspl_twsit, cspl_twist)

    chord = cspl.(rbe)
    twist = tspl.(rbe) 

    sweep = zeros(nrow)
    dihed = zeros(nrow)
    n_seg = repeat([nseg],nrow)
    type_seg = repeat(["uniform"],nrow)
    airfoil = repeat([afname],nrow)
    geom = hcat(collect(rbe),chord, twist,sweep,dihed,n_seg,type_seg,airfoil)

    return geom 
end 

#------直升机相关简单数学公式------
# 悬停效率计算公式
hoverFM(R,T,P;rho=1.225) = T*sqrt(T/(2*rho*pi*R^2))/P 
hoverFM(R,T,Q,Ω;rho=1.225) = T*sqrt(T/(2*rho*pi*R^2))/(Q*Ω)

# 挥舞铰惯性矩
function moi_flap(dm, rr, ef)
    return 1.0/3.0*dm*rr^3*(1-ef)^3
end 
  
# 挥舞铰质量静矩
function ms_flap(dm, rr, ef)
    return 0.5*dm*rr^2*(1-ef)^2
end 

# 桨叶挥舞洛克数 γ 经验公式
function lock_num(cla, chord, rr, dm, ef;rho=1.125)
    return rho*cla*chord*rr^4/moi_flap(dm, rr, ef)
end

# Prouty的桨尖损失经验公式
function TipLoss(ct; nb=4) #[Prouty] 

    if ct <= 0
        tiploss = 1.0
    elseif ct <= 0.006
        tiploss = 1 - 0.06/nb
    else 
        tiploss = 1 - sqrt(4.27*ct-0.01)/nb 
    end 

    return tiploss 
end 

end # module
