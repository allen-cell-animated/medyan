function dist!(out,x)
	d_low = 6.0e3
	d_high = 7.0e3
	N = size(x,2)
	counter = 1
	allpairs = 1
	# C = zeros(Float64,N*(N-1)/2 |> Int)
	for i in 1:N
		for j in i+1:N
			 dx = x[1,i]-x[1,j]
			 dy = x[2,i]-x[2,j]
			 dz = x[3,i]-x[3,j]
			c = dx*dx+dy*dy+dz*dz
			# @show i, j, c
			# println(dx," ", dy, " ", dz)
			# C[allpairs] = c
			if d_low < c < d_high
				out[1,counter] = i
				out[2,counter] =j
				counter+=1
			end
			allpairs+=1
		end
	end
	# return C
	return allpairs-1, counter-1
end

function findpairs()
	N = 1500
	x=100*rand(3,N) 
	out_N = N*(N-1)/2 |> Int
	out = zeros(Int,2,out_N)
	# C = dist!(out,x)
	allpairs, found_pairs = @time dist!(out,x)
	@show allpairs, found_pairs
	return found_pairs, out
end

found_pairs, out = findpairs()
found_pairs, out = findpairs()


