module RoTools

import Interpolations
import BSplines

export BladeGeom, itr2d, hoverFM

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

function BladeGeom(mat::Matrix)
    r = mat[:,1] # radius
    c = mat[:,2] # chord
    t = mat[:,3] # twist 

    itr_c = itr2d(r, c)
    itr_t = itr2d(r, t)

    rbe = range(r[1],r[end],20)
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

# 悬停效率计算公式
hoverFM(R,T,P;rho=1.225) = T*sqrt(T/(2*rho*pi*R^2))/P 
hoverFM(R,T,Q,Ω;rho=1.225) = T*sqrt(T/(2*rho*pi*R^2))/(Q*Ω)

end # module
