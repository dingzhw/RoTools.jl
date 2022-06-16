module RoTools

using Interpolations

export BladeGeom 

# 桨叶弦长和扭转分布规范化
function itr2d(x,y; extrapolation_method = Line()) # make convience LinearInterpolation of 2d array 
    res = LinearInterpolation(x,y,extrapolation_bc=extrapolation_method)
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

end # module
