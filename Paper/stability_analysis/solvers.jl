# Complementary code to paper TODO
# Contains the code for the solvers used in examples.jl
module Solvers

import LinearAlgebra
const LA=LinearAlgebra

struct Data

    currintensity::Vector{Float64}
    allintensities::Array{Float64}
    backgroundintensity::Array{Float64}#for the boundary conditions
    χ::Array{Float64}
    η::Array{Float64}
    npoints::Int
    nfreqs::Int
    v::Vector{Float64}#implied to be divided by c the speed of light
    x::Vector{Float64}
    ν::Array{Float64}
    S::Float64#For these tests, the source function S will always be constant

end

function data(backgroundintensity, x, v, χ, η, ν, S)
    @views currintensity=copy(backgroundintensity[:,1])
    nfreqs,_=size(backgroundintensity)
    npoints=length(x)
    allintensities=zeros(Float64, nfreqs, npoints)
    @views allintensities[:,1].=currintensity#backgroundintensity[:,1]

    return Data(currintensity, allintensities, backgroundintensity, χ, η, npoints, nfreqs, v, x, ν, S)
end


#computes intensity for a single distance increment using first order fully explicit discretization
function firstorderexplicit(previndex, data)
    Δx=data.x[previndex+1]-data.x[previndex]#assumes x strictly increasing
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ

    currintensity[1:nfreqs].+=Δx.*(η[1:nfreqs, previndex]-currintensity[1:nfreqs].*χ[1:nfreqs, previndex])

    data.allintensities[:,previndex+1]=currintensity;

    return
end

function firstorderexplicitsolver(data::Data)
    for i ∈ 1:data.npoints-1
        firstorderexplicit(i, data)
    end

    return
end


#computes intensity for a single optical depth increment using second order semi implicit discretization (also second order for the frequency derivative)
# Also adds approximate inner boundary conditions if the frequency difference were to change sign due to changing line width
function secondordershortcharmoving(previndex, data)
    Δx=(data.x[previndex+1]-data.x[previndex])#assumes x strictly increasing
    Δxdiv2=(data.x[previndex+1]-data.x[previndex])/2.0#assumes x strictly increasing
    Δvdiv2=(data.v[previndex+1]-data.v[previndex])/2.0

    @views Δνdiv2=(data.ν[:,previndex+1]-data.ν[:,previndex])/2.0
    @views Δν=(data.ν[:,previndex+1]-data.ν[:,previndex])
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ
    ν=data.ν
    S=η./χ
    # S=fill(data.S, size(η))
    #bounding Δτ from below, as division by almost 0 errors can happen. (we will be dividing by Δτ.^2)
    Δτ=Δx*(χ[:,previndex+1].+χ[:,previndex])/2.0
    minτ=1.0E-10
    toolow = findall(Δτ.<=minτ)
    Δτ[toolow] .= minτ
    Δτdiv2 = Δτ./2.0

    expminτdiv2=exp.(-Δτdiv2)
    expminτ=exp.(-2.0.*Δτdiv2)
    onemexpminτdiv2=-expm1.(-Δτdiv2)
    onemexpminτ=-expm1.(-2.0 .*Δτdiv2)

    # lineν=data.lineν[previndex]
    # nextlineν=data.lineν[previndex+1]

    #now compute for which frequency points we need a forward or backward discretization
    currν=data.ν[:, previndex]
    nextν=data.ν[:, previndex+1]

    # forwardfreqdisc=(nextν.-currν.+nextlineν.-lineν.>0.0)
    forwardfreqdisc=(nextν.-currν.>0.0)


    starting_upwind=forwardfreqdisc[1]
    ending_upwind=forwardfreqdisc[nfreqs]

    other_discretization_direction=[!(starting_upwind==forwardfreqdisc[index]) for index in 1:length(forwardfreqdisc)]

    inflection_point=findfirst(other_discretization_direction)

    #Default discretization direction
    if isnothing(inflection_point)
        if starting_upwind
            boundary_points=[(nfreqs-1,nfreqs)]#is tuple, as they should be treated together
            boundary_point_is_outer=[true]
            upwind_points=1:(nfreqs-2)
            downwind_points=[]
        else
            boundary_points=[(1,2)]
            boundary_point_is_outer=[true]
            upwind_points=[]
            downwind_points=3:nfreqs
        end
    else
        if starting_upwind
            boundary_points=[(inflection_point-1,inflection_point)]
            boundary_point_is_outer=[false]
            #thus only inner boundary points
            upwind_points=1:(inflection_point-2)
            downwind_points=inflection_point+1:nfreqs
        else
            boundary_points=[(1,2),(nfreqs-1,nfreqs)]
            boundary_point_is_outer=[true, true]
            #thus only outer boundary points
            upwind_points=inflection_point:nfreqs-2
            downwind_points=3:(inflection_point-1)
        end
    end

    #forward discretization

    Δνsmall=(view(ν, upwind_points.+2, previndex)-view(ν, upwind_points.+1, previndex))
    Δνlarge=(view(ν, upwind_points.+2, previndex)-view(ν, upwind_points, previndex))
    a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
    b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
    c=.-a.-b;

    #is the numerical derivative of I
    Δνterm=(a.*view(currintensity, upwind_points.+2).+b.*view(currintensity, upwind_points.+1).+c.*view(currintensity, upwind_points))#.*lineν

    @views sourceterm=((onemexpminτ[upwind_points].-Δτ[upwind_points])./Δτ[upwind_points].+onemexpminτ[upwind_points]).*S[upwind_points, previndex].+(Δτ[upwind_points].-onemexpminτ[upwind_points]).*S[upwind_points, previndex+1]./Δτ[upwind_points]
    @views curr_factor=onemexpminτ[upwind_points]./Δτ[upwind_points].+(onemexpminτ[upwind_points].-Δτ[upwind_points])./Δτ[upwind_points].^2

    @views currintensity[upwind_points]=currintensity[upwind_points].*expminτ[upwind_points].+sourceterm.+curr_factor.*Δνterm.*Δν[upwind_points]

    #end forward discretization explicit part
    #now do backward discretization explicit part

    Δνsmall=(view(ν, downwind_points.-1, previndex)-view(ν, downwind_points.-2, previndex))
    Δνlarge=(view(ν, downwind_points, previndex)-view(ν, downwind_points.-2, previndex))
    a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
    b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
    c=.-a.-b;

    Δνterm=(a.*view(currintensity, downwind_points.-2).+b.*view(currintensity, downwind_points.-1).+c.*view(currintensity, downwind_points))

    @views sourceterm=((onemexpminτ[downwind_points].-Δτ[downwind_points])./Δτ[downwind_points].+onemexpminτ[downwind_points]).*S[downwind_points, previndex].+(Δτ[downwind_points].-onemexpminτ[downwind_points]).*S[downwind_points, previndex+1]./Δτ[downwind_points]
    @views curr_factor=onemexpminτ[downwind_points]./Δτ[downwind_points].+(onemexpminτ[downwind_points].-Δτ[downwind_points])./Δτ[downwind_points].^2

    @views currintensity[downwind_points]=currintensity[downwind_points].*expminτ[downwind_points].+sourceterm.-curr_factor.*Δνterm.*Δν[downwind_points]

    #end backward discretization explicit part

    #now apply approximate inner boundary conditions by interpolating and applying the formal solution
    for bdy_tpl_index ∈ 1:length(boundary_points)
        indices=[boundary_points[bdy_tpl_index][i] for i ∈ 1:length(boundary_points[bdy_tpl_index])]
        if (boundary_point_is_outer[bdy_tpl_index])
            #then just set it to the boundary value
            @views currintensity[indices]=bdyintensity[indices, previndex+1]
        else
            νdiff=currν[indices[2]]-currν[indices[1]]
            νleft=currν[indices[1]]
            #interpolate them to their next frequencies (linearly)
            @views currintensity[indices]=currintensity[indices[1]].+
                    (nextν[indices].-νleft)./νdiff.*(currintensity[indices[2]].-currintensity[indices[1]])
            #also interpolate some opacities and emissivities
            @views ηint=η[indices[1], previndex].+(nextν[indices].-νleft)./νdiff.*(η[indices[2], previndex].-η[indices[1], previndex])
            @views χint=χ[indices[1], previndex].+(nextν[indices].-νleft)./νdiff.*(χ[indices[2], previndex].-χ[indices[1], previndex])
            #now apply 2nd order static solver to them
            @views currintensity[indices]=((currintensity[indices].*(1.0 .-Δxdiv2.*χint).+Δxdiv2.*(ηint.+η[indices, previndex+1]))
                                           ./(1.0.+Δxdiv2.*χ[indices, previndex+1]))
        end
    end

    # lineν=data.lineν[previndex+1]

    #upwind implicit part
    #err, just ignore implicit part for now if no points need to be computed with this discretization
    if length(upwind_points)>0

        Δνsmall=(view(ν, upwind_points.+2, previndex+1)-view(ν, upwind_points.+1, previndex+1))
        Δνlarge=(view(ν, upwind_points.+2, previndex+1)-view(ν, upwind_points, previndex+1))
        a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
        b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
        c=.-a.-b;

        Δνnext=(view(ν, upwind_points, previndex+1)-view(ν, upwind_points, previndex))

        @views curr_factor=(Δτ[upwind_points].-onemexpminτ[upwind_points])./Δτ[upwind_points].^2

        matrixsize=length(upwind_points)+2#+2 boundary conditions at the end
        #setting up the matrix
        diagonal=ones(matrixsize)
        offdiagonal=zeros(matrixsize-1)
        secondoffdiagonal=zeros(matrixsize-2)

        @views diagonal[1:(matrixsize-2)].+=-curr_factor.*Δν[upwind_points].*c
        @views offdiagonal[1:(matrixsize-2)]=-curr_factor.*Δν[upwind_points].*b
        @views secondoffdiagonal=-curr_factor.*Δν[upwind_points].*a

        #inefficient, as julia stores the entire matrix, but this should work
        matrix=LA.diagm(0 => diagonal,1=>offdiagonal, 2=>secondoffdiagonal)

        rangeincludingbdy=UnitRange(first(upwind_points), last(upwind_points)+2)
        @views currintensity[rangeincludingbdy].=(matrix \ currintensity[rangeincludingbdy])

    end

    #end upwind implicit part

    #start downwind implicit part
    #err, just ignore implicit part for now if no points need to be computed with this discretization
    if length(downwind_points)>0

        Δνsmall=(view(ν, downwind_points.-1, previndex+1)-view(ν, downwind_points.-2, previndex+1))
        Δνlarge=(view(ν, downwind_points, previndex+1)-view(ν, downwind_points.-2, previndex+1))
        a=-Δνsmall./(Δνlarge.^2 .-Δνsmall.*Δνlarge)
        b=Δνlarge./(Δνsmall.*Δνlarge.-Δνsmall.^2)
        c=.-a.-b;

        Δνnext=(view(ν, downwind_points, previndex+1)-view(ν, downwind_points, previndex))

        @views curr_factor=(Δτ[downwind_points].-onemexpminτ[downwind_points])./Δτ[downwind_points].^2

        matrixsize=length(downwind_points)+2#+2 boundary conditions at the beginning
        #setting up the matrix
        diagonal=ones(matrixsize)
        offdiagonal=zeros(matrixsize-1)
        secondoffdiagonal=zeros(matrixsize-2)

        #note: FD formula for backward is exactly the same, except we need to change sign of coefficients

        @views diagonal[3:matrixsize].+=curr_factor.*Δν[downwind_points].*c
        @views offdiagonal[2:(matrixsize-1)]=curr_factor.*Δν[downwind_points].*b
        @views secondoffdiagonal=curr_factor.*Δν[downwind_points].*a

        #inefficient, as julia stores the entire matrix, but this should work
        matrix=LA.diagm(0 => diagonal,-1=>offdiagonal, -2=>secondoffdiagonal)
        rangeincludingbdy=UnitRange(first(downwind_points)-2, last(downwind_points))

        @views currintensity[rangeincludingbdy].=(matrix \ currintensity[rangeincludingbdy])

    end

    data.allintensities[:,previndex+1]=currintensity;
    return
end


function secondordershortcharmovingsolver(data::Data)
    for i ∈ 1:data.npoints-1
        secondordershortcharmoving(i, data)
    end

    return
end


#computes intensity using second order semi implicit discretization (also second order for the frequency derivative). Matches the frequency indices to obtain better stability.
function secondordermovingshortcharfreqmatch(previndex, data)
    Δx=(data.x[previndex+1]-data.x[previndex])#assumes x strictly increasing
    Δxdiv2=(data.x[previndex+1]-data.x[previndex])/2.0#assumes x strictly increasing
    Δvdiv2=(data.v[previndex+1]-data.v[previndex])/2.0
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ
    ν=data.ν
    S=η./χ

    #now compute for which frequency points we need a forward or backward discretization
    currν=data.ν[:, previndex]
    nextν=data.ν[:, previndex+1]

    forwardfreqdisc=(nextν.-currν.>0.0)
    #we can impose the direction of discretizing for this solver, this is done by majority vote
    forwarddisc=(count(x->x, forwardfreqdisc)>=nfreqs/2)

    #for all points, determine the interval of the previous freqs
    #assume frequencies sorted
    nextfreqidx=1
    currfreqidx=1
    is_boundary_point=trues(nfreqs)#BitArray(true, nfreqs)
    curr_corr_idx=zeros(Int, nfreqs)#contains the corresponding closest frequency indices (at previndex) to each frequency at previndex+1
    #add extra boundary conditions due to mismatching
    while ((currfreqidx <= nfreqs) && (nextfreqidx <= nfreqs))
        if (currν[currfreqidx]<nextν[nextfreqidx])
            currfreqidx+=1
        else
            is_boundary_point[nextfreqidx]=false
            if forwarddisc
                curr_corr_idx[nextfreqidx]=currfreqidx-1
            else
                curr_corr_idx[nextfreqidx]=currfreqidx
            end
            nextfreqidx+=1
        end
    end
    #and set boundary condition if left boundary point would be 0 (i.e. outside domain)
    for tempidx∈1:nfreqs
        if curr_corr_idx[tempidx]==0
            is_boundary_point[tempidx]=true
        end
    end
    #also set the outmost point(s) in the direction to be (a) boundary point(s).
    if forwarddisc
        is_boundary_point[nfreqs-1]=true
        is_boundary_point[nfreqs]=true
        is_boundary_point[curr_corr_idx.==nfreqs-1].=true
        is_boundary_point[curr_corr_idx.==nfreqs].=true
    else
        is_boundary_point[1]=true
        is_boundary_point[2]=true
        is_boundary_point[curr_corr_idx.==1].=true
        is_boundary_point[curr_corr_idx.==2].=true
    end

    #forward discretization
    boundary_points=findall(is_boundary_point.==true)
    if forwarddisc
        upwind_points=findall(is_boundary_point.==false)
        upwind_curr_corr_idx=curr_corr_idx[upwind_points]

        Δτ=Δx*(χ[upwind_points,previndex+1].+χ[upwind_curr_corr_idx,previndex])/2.0
        minτ=1.0E-10
        toolow = findall(Δτ.<=minτ)
        Δτ[toolow] .= minτ
        Δτdiv2 = Δτ./2.0

        expminτdiv2=exp.(-Δτdiv2)
        expminτ=exp.(-2.0.*Δτdiv2)
        onemexpminτdiv2=-expm1.(-Δτdiv2)
        onemexpminτ=-expm1.(-2.0.*Δτdiv2)

        @views Δνdiv2=(data.ν[upwind_points,previndex+1]-data.ν[upwind_curr_corr_idx,previndex])/2.0
        @views Δν=(data.ν[upwind_points,previndex+1]-data.ν[upwind_curr_corr_idx,previndex])

        Δνsmall=(view(ν, upwind_curr_corr_idx.+2, previndex)-view(ν, upwind_curr_corr_idx.+1, previndex))
        Δνlarge=(view(ν, upwind_curr_corr_idx.+2, previndex)-view(ν, upwind_curr_corr_idx, previndex))
        a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
        b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
        c=.-a.-b;

        #is the numerical derivative of I
        Δνterm=(a.*view(currintensity, upwind_curr_corr_idx.+2).+b.*view(currintensity, upwind_curr_corr_idx.+1).+c.*view(currintensity, upwind_curr_corr_idx))#.*lineν

        @views sourceterm=((onemexpminτ.-Δτ)./Δτ.+onemexpminτ).*S[upwind_curr_corr_idx, previndex].+(Δτ.-onemexpminτ).*S[upwind_points, previndex+1]./Δτ
        @views curr_factor=onemexpminτ./Δτ.+(onemexpminτ.-Δτ)./Δτ.^2
        @views currintensity[upwind_points]=currintensity[upwind_curr_corr_idx].*expminτ.+sourceterm.+curr_factor.*Δνterm.*Δν

        #end forward discretization explicit part
    else
        #do backward discretization explicit part
        downwind_points=findall(is_boundary_point.==false)

        downwind_curr_corr_idx=curr_corr_idx[downwind_points]

        Δτ=Δx*(χ[downwind_points,previndex+1].+χ[downwind_curr_corr_idx,previndex])/2.0
        minτ=1.0E-10
        toolow = findall(Δτ.<=minτ)
        Δτ[toolow] .= minτ
        Δτdiv2 = Δτ./2.0

        expminτdiv2=exp.(-Δτdiv2)
        expminτ=exp.(-2.0.*Δτdiv2)
        onemexpminτdiv2=-expm1.(-Δτdiv2)
        onemexpminτ=-expm1.(-2.0.*Δτdiv2)

        @views Δνdiv2=(data.ν[downwind_points,previndex+1]-data.ν[downwind_curr_corr_idx,previndex])/2.0
        @views Δν=(data.ν[downwind_points,previndex+1]-data.ν[downwind_curr_corr_idx,previndex])

        Δνsmall=(view(ν, downwind_curr_corr_idx.-1, previndex)-view(ν, downwind_curr_corr_idx.-2, previndex))
        Δνlarge=(view(ν, downwind_curr_corr_idx, previndex)-view(ν, downwind_curr_corr_idx.-2, previndex))
        a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
        b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
        c=.-a.-b;

        #is the numerical derivative of I
        Δνterm=(a.*view(currintensity, downwind_curr_corr_idx.-2).+b.*view(currintensity, downwind_curr_corr_idx.-1).+c.*view(currintensity, downwind_curr_corr_idx))#.*lineν

        @views sourceterm=((onemexpminτ.-Δτ)./Δτ.+onemexpminτ).*S[downwind_curr_corr_idx, previndex].+(Δτ.-onemexpminτ).*S[downwind_points, previndex+1]./Δτ
        @views curr_factor=onemexpminτ./Δτ.+(onemexpminτ.-Δτ)./Δτ.^2
        @views currintensity[downwind_points]=currintensity[downwind_curr_corr_idx].*expminτ.+sourceterm.+curr_factor.*Δνterm.*-Δν

        #end backward discretization explicit part
    end

    #end backward discretization explicit part

    #now apply boundary conditions
    @views currintensity[boundary_points]=bdyintensity[boundary_points, previndex+1]

    #upwind implicit part
    if forwarddisc

        upwind_points=findall(is_boundary_point.==false)
        upwind_curr_corr_idx=curr_corr_idx[upwind_points]

        Δνsmall=(view(ν, upwind_points.+2, previndex+1)-view(ν, upwind_points.+1, previndex+1))
        Δνlarge=(view(ν, upwind_points.+2, previndex+1)-view(ν, upwind_points, previndex+1))
        a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
        b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
        c=.-a.-b;

        Δτ=Δx*(χ[upwind_points,previndex+1].+χ[upwind_curr_corr_idx,previndex])/2.0
        minτ=1.0E-10
        toolow = findall(Δτ.<=minτ)
        Δτ[toolow] .= minτ
        Δτdiv2 = Δτ./2.0

        expminτdiv2=exp.(-Δτdiv2)
        expminτ=exp.(-2.0.*Δτdiv2)
        onemexpminτdiv2=-expm1.(-Δτdiv2)
        onemexpminτ=-expm1.(-2.0.*Δτdiv2)

        Δν=(view(ν, upwind_points, previndex+1)-view(ν, upwind_curr_corr_idx, previndex))

        @views curr_factor=(Δτ.-onemexpminτ)./Δτ.^2

        matrixsize=length(upwind_points)+2#+2 boundary conditions at the end
        #setting up the matrix
        diagonal=ones(matrixsize)
        offdiagonal=zeros(matrixsize-1)
        secondoffdiagonal=zeros(matrixsize-2)

        @views diagonal[1:(matrixsize-2)].+=-curr_factor.*Δν.*c
        @views offdiagonal[1:(matrixsize-2)]=-curr_factor.*Δν.*b
        @views secondoffdiagonal=-curr_factor.*Δν.*a

        #inefficient, as julia stores the entire matrix, but this should work
        matrix=LA.diagm(0 => diagonal,1=>offdiagonal, 2=>secondoffdiagonal)

        rangeincludingbdy=UnitRange(first(upwind_points), last(upwind_points)+2)
        @views currintensity[rangeincludingbdy].=(matrix \ currintensity[rangeincludingbdy])

        #end upwind implicit part
    else

        #start downwind implicit part

        downwind_points=findall(is_boundary_point.==false)

        downwind_curr_corr_idx=curr_corr_idx[downwind_points]

        Δτ=Δx*(χ[downwind_points,previndex+1].+χ[downwind_curr_corr_idx,previndex])/2.0
        minτ=1.0E-10
        toolow = findall(Δτ.<=minτ)
        Δτ[toolow] .= minτ
        Δτdiv2 = Δτ./2.0

        expminτdiv2=exp.(-Δτdiv2)
        expminτ=exp.(-2.0.*Δτdiv2)
        onemexpminτdiv2=-expm1.(-Δτdiv2)
        onemexpminτ=-expm1.(-2.0.*Δτdiv2)

        Δν=(view(ν, downwind_points, previndex+1)-view(ν, downwind_curr_corr_idx, previndex))

        Δνsmall=(view(ν, downwind_points.-1, previndex+1)-view(ν, downwind_points.-2, previndex+1))
        Δνlarge=(view(ν, downwind_points, previndex+1)-view(ν, downwind_points.-2, previndex+1))
        a=-Δνsmall./(Δνlarge.^2 .-Δνsmall.*Δνlarge)
        b=Δνlarge./(Δνsmall.*Δνlarge.-Δνsmall.^2)
        c=.-a.-b;

        @views curr_factor=(Δτ.-onemexpminτ)./Δτ.^2

        matrixsize=length(downwind_points)+2#+2 boundary conditions at the beginning
        #setting up the matrix
        diagonal=ones(matrixsize)
        offdiagonal=zeros(matrixsize-1)
        secondoffdiagonal=zeros(matrixsize-2)

        #note: FD formula for backward is exactly the same, except we need to change sign of coefficients

        println("factor on diagonal: ", curr_factor.*Δν.*c)

        @views diagonal[3:matrixsize].+=curr_factor.*Δν.*c
        @views offdiagonal[2:(matrixsize-1)]=curr_factor.*Δν.*b
        @views secondoffdiagonal=curr_factor.*Δν.*a

        #inefficient, as julia stores the entire matrix, but this should work
        matrix=LA.diagm(0 => diagonal,-1=>offdiagonal, -2=>secondoffdiagonal)

        rangeincludingbdy=UnitRange(first(downwind_points)-2, last(downwind_points))

        @views currintensity[rangeincludingbdy].=(matrix \ currintensity[rangeincludingbdy])

    end

    data.allintensities[:,previndex+1]=currintensity;
    return
end


function secondordermovingshortcharfreqmatchsolver(data::Data)
    for i ∈ 1:data.npoints-1
        secondordermovingshortcharfreqmatch(i, data)
    end

    return
end


# second order feautrier solver for comparison
# computes mean intensities u on an entire ray for a single frequency at a time
# Uses notation from rybicki & hummer 1991
function second_order_Feautrier(data::Data)

    if data.nfreqs≠1
        error("nfreqs!=1")
    end

    N = data.npoints
    #everywhere, the source is constant in these examples
    S = ones(N).*data.S
    dtau = (data.χ[1,1:N-1]+data.χ[1,2:N]) .* (data.x[2:N]-data.x[1:N-1]) ./2.0#defined in range 1:N-1
    I0 = data.backgroundintensity[1,1]
    In = data.backgroundintensity[1,N]

    A = zeros(N)#practically from 2:N
    B = zeros(N)
    C = zeros(N)#practically from 1:N-1
    F = zeros(N)
    z = zeros(N)
    u = zeros(N)

    A[2:N-1] = 2.0 ./ ((dtau[1:N-2] .+ dtau[2:N-1]) .* dtau[1:N-2])
    C[2:N-1] = 2.0 ./ ((dtau[1:N-2] .+ dtau[2:N-1]) .* dtau[2:N-1])

    B = 1.0 .+ A .+ C

    B[1] = 1.0 + 2.0/dtau[1] + 2.0/dtau[1]^2
    C[1] =                     2.0/dtau[1]^2

    S[1] += 2.0/dtau[1] * I0

    B[N] = 1.0 + 2.0/dtau[N-1] + 2.0/dtau[N-1]^2
    A[N] =                       2.0/dtau[N-1]^2

    Bn_min_An = 1.0 + 2.0/dtau[N-1]

    S[N] += 2.0/dtau[N-1] * In

    F[1] = B[1]/C[1] - 1.0
    z[1] = S[1]/B[1]

    # forward elimination
    for n ∈ 2:N
        F[n] = (1.0 + A[n]*F[n-1]/(1.0 + F[n-1])) / C[n]
        z[n] = (S[n] + A[n]*z[n-1])/(C[n]*(1.0 + F[n]))
    end

    u[N] = (S[N] + A[N] * z[N-1]) / (B[N] * F[N-1] + Bn_min_An) * (1.0 + F[N-1])
    data.allintensities[1,N]=u[N]

    # backward elimination
    for n ∈ N-1:-1:1
        u[n] = z[n] + u[n+1] / (1.0 + F[n])
        data.allintensities[1,n]=u[n]
    end

end


# fourth order feautrier solver
# computes mean intensities u on an entire ray for a single frequency at a time
# Uses notation from rybicki & hummer 1991
function fourth_order_Feautrier(data::Data)

    if data.nfreqs≠1
        error("nfreqs!=1")
    end

    N = data.npoints
    #everywhere, the source is constant in these examples
    #because the source is constant everywhere, we will not bother correctly summing it with the nearby source functions with the correct weights (as they sum to 1)
    S = ones(N).*data.S
    dtau = (data.χ[1,1:N-1]+data.χ[1,2:N]) .* (data.x[2:N]-data.x[1:N-1]) ./2.0#defined in range 1:N-1
    I0 = data.backgroundintensity[1,1]
    In = data.backgroundintensity[1,N]

    A = zeros(N)#practically from 2:N
    B = zeros(N)
    C = zeros(N)#practically from 1:N-1
    F = zeros(N)
    z = zeros(N)
    u = zeros(N)

    # default feautrier coefficients
    A[2:N-1] = 2.0 ./ ((dtau[1:N-2] .+ dtau[2:N-1]) .* dtau[1:N-2])
    C[2:N-1] = 2.0 ./ ((dtau[1:N-2] .+ dtau[2:N-1]) .* dtau[2:N-1])

    #slight change due to fourth order
    A[2:N-1] = A[2:N-1] .- 1.0./12.0 .* (2.0 .- dtau[2:N-1].^2 .* A[2:N-1])
    C[2:N-1] = C[2:N-1] .- 1.0./12.0 .* (2.0 .- dtau[1:N-2].^2 .* C[2:N-1])

    B = 1.0 .+ A .+ C

    B[1] = 2.0/3.0 + 2.0/dtau[1] + 2.0/dtau[1]^2
    C[1] = -1.0/3.0              + 2.0/dtau[1]^2

    S[1] += 2.0/dtau[1] * I0

    B[N] = 2.0/3.0 + 2.0/dtau[N-1] + 2.0/dtau[N-1]^2
    A[N] = -1.0/3.0                + 2.0/dtau[N-1]^2

    Bn_min_An = 1.0 + 2.0/dtau[N-1]

    S[N] += 2.0/dtau[N-1] * In

    F[1] = B[1]/C[1] - 1.0
    z[1] = S[1]/B[1]

    # forward elimination
    for n ∈ 2:N
        F[n] = (1.0 + A[n]*F[n-1]/(1.0 + F[n-1])) / C[n]
        z[n] = (S[n] + A[n]*z[n-1])/(C[n]*(1.0 + F[n]))
    end

    u[N] = (S[N] + A[N] * z[N-1]) / (B[N] * F[N-1] + Bn_min_An) * (1.0 + F[N-1])
    data.allintensities[1,N]=u[N]

    # backward elimination
    for n ∈ N-1:-1:1
        u[n] = z[n] + u[n+1] / (1.0 + F[n])
        data.allintensities[1,n]=u[n]
    end

end


end
