#=
    Given a function that is piecewise constant, and changes 
    values at specified times, evaluate the function at the newTimes.

Parameters
----------
curtimes: array
         Times that the values of the trajectory change at. Assumes starts at t=0.

traject: array
         Values of the trajectory after each time change in curtimes. Assumes values >= 0.

newtimes: array
         Times to evaluate the trajectory at.

Returns
-------
newtraject: array
         Value of the trajectory at the times in newtimes
=#
function get_value_at(curtimes, traject, newtimes)
    newtraject = zeros( length(newtimes) )
    idx        = 1

    for i in eachindex(curtimes)[2:end]
        while newtimes[idx] < curtimes[i] # If newtimes[idx] > curtimes[i], then skip this i and move on to next i in the for loop
            newtraject[idx] = traject[i-1]
            idx = idx + 1
            if idx > length(newtimes)
                return newtraject
            end
        end
    end

    # at this point, we have go through curtimes. However, if some of the new times are > the last original time,
    # just take the last value from curTimes
    if idx <= length(newtimes)
        newtraject[idx:end] .= traject[end]
    end

    return newtraject
end
