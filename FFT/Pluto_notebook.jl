### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# ╔═╡ 85babbf4-3aab-11eb-1418-2140dcd88435

using DelimitedFiles	


# ╔═╡ 0b1acabe-3ab6-11eb-3141-371834a04fbd
dir_str = readdir()

# ╔═╡ 97a4c138-3ab6-11eb-0365-d7c2e6a7ffe5
dir_str[1]

# ╔═╡ 50e7b37a-3ab6-11eb-184b-71e50ff72c3d
ind_Hz = findall(contains.(dir_str, "Hz"))

# ╔═╡ afb55dae-3ab7-11eb-04fb-31cfab44f0bc
function substr_extract(original_string, begin_string, end_string)
	    ind_desc = findfirst(begin_string, original_string); 
        name_str = ind_desc[end]+1;
		ind_desc1 = findnext(end_string, original_string, name_str);
        name_end = ind_desc1[1]-1;
        Hz = original_string[name_str:name_end];		
return parse(Float64, Hz)
end

# ╔═╡ adac1c30-3ac3-11eb-2d7a-c34e1278251f
function substr_extract1(original_string, begin_string, end_string)
	    ind_desc = findfirst(begin_string, original_string); 
        name_str = ind_desc[end]+1;
		ind_desc1 = findnext(end_string, original_string, name_str);
        name_end = ind_desc1[1]-1;
        Hz = original_string[name_str:name_end];		
return Hz
end

# ╔═╡ 5cf4d328-3ab8-11eb-1948-5faeb828eca8
length(ind_Hz)

# ╔═╡ 071e8caa-3ab8-11eb-1e4f-b159d87cabe1
#begin
#for i=1:length(ind_Hz)
substr_extract(dir_str[ind_Hz[38]], "freq","Hz")
#end
#end

# ╔═╡ 50830bfc-3ab6-11eb-0c58-bfa15d7747d0
Hz_data = [substr_extract(dir_str[ind_Hz[i]],"freq", "Hz") for i =1:length(ind_Hz)]


# ╔═╡ 5021084e-3ab6-11eb-15d3-279d2c4b11f0
Hz_data1 = sort(Hz_data)

# ╔═╡ 32cd8376-3ac0-11eb-2061-0f28dc1a9e77
ind_Int = findall(iszero, Hz_data1 -floor.(Hz_data1))

# ╔═╡ c8666f5a-3ac1-11eb-0025-95a05abc171d
typeof(ind_Int)

# ╔═╡ 57c399e4-3ac1-11eb-3249-91634004453c
ind_Float = findall(!iszero, Hz_data1-floor.(Hz_data1))

# ╔═╡ 7f4d6a52-3aa9-11eb-06c9-c14eb7146d4b

file_folders = ["./freq"*string(Int(Hz_data1[i]))*"Hz/t3s2_type4_inj.inf" for i in ind_Int]
#
	


# ╔═╡ df33fbb2-3ac1-11eb-1fe7-2f289db520b5
file_folders2= [ "./freq"*string(Hz_data1[i])*"Hz/t3s2_type4_inj.inf" for i in ind_Float]

# ╔═╡ 3e88fd4c-3ac2-11eb-354d-f19f982f34aa
append!(file_folders, file_folders2)

# ╔═╡ 49a76a08-3ac2-11eb-0043-d7fba49a40ea
length(file_folders)

# ╔═╡ 1a590fa6-3ac3-11eb-388f-1d68e504ad14
file_folders[1]

# ╔═╡ 81c43378-3ac3-11eb-02c1-8b3f25d8b640
begin
file_folder = file_folders[1]
ff = open(file_folder)
data_ff =readlines(ff)
names1 =[substr_extract1(data_ff[i], "Desc=\"","\"") for i in 1:length(data_ff)]
end

# ╔═╡ 505751d2-3aab-11eb-1642-c55919d3ea8f
names1

# ╔═╡ 6aa2a052-3aab-11eb-2420-7bb2200dfb25
N_files = Int(floor(length(names1)/10)+1)

# ╔═╡ 6c54fca0-3ac5-11eb-3c83-31e7a826392f
Hz_string = [string(Hz_data1[i]) for i in ind_Float]

# ╔═╡ c0d9925c-3ac5-11eb-2a23-fd1ab847886b
Hz_string1 = [string(Int(Hz_data1[i])) for i in ind_Int]

# ╔═╡ bab146bc-3ac5-11eb-2b8d-ff2fe74d6ba3
append!(Hz_string, Hz_string1)

# ╔═╡ 7502cebe-3aab-11eb-1318-e5aea2a6c9da
filenames = ["./freq"*Hz_string[j]*"Hz/t3s2_type4_inj_"*string(i, pad=2)*".out" for i=1:N_files, j=1:length(Hz_string)]

# ╔═╡ 37dce552-3ac6-11eb-19e0-cd99abf02a5c
size(filenames)

# ╔═╡ 9b61420a-3ac6-11eb-0fd4-6f845ad884ae
k = 5; #  freq index

# ╔═╡ ac95f44e-3ac6-11eb-3a09-b77cb8f5e850
Hz_string[k]

# ╔═╡ 4727a398-3ab0-11eb-35b6-5d4d0b4fb102
data_many =readdlm.(filenames[:,k],  header=false)

# ╔═╡ 8e3d5994-3aab-11eb-1daa-77db968d33c4
begin
function combine_datafiles(data_many1, n)
	xx=data_many1[1];
	for i in 2:n
		xx=hcat(xx, data_many1[i][:,2:end]);
	end
	return xx;
end
Data_total = combine_datafiles(data_many,N_files);
end

# ╔═╡ fe1ca76a-3aa8-11eb-1259-a1d5fdb20f4a
begin
	
using Plots
#using GKS
plot(Data_total[:,1], Data_total[:,2:20], title = Hz_string[k]*" Hz perturbation")

# time interval is 0.0005 seconds or 2000 Hz

using FFTW
using AbstractFFTs

tt = Data_total[:,1];
Ts = tt[2]-tt[1];
# Start time
t0 = 1.1; 
N = 7800-1;
N_start = Int(floor(t0/Ts));
tmax = t0 + N * Ts
N_end = N_start + N;
	
nn_start = Int(floor(0.9/Ts));

# time coordinate
#t = [t0:Ts:tmax;]
t = Data_total[N_start:N_end,1]
# signal
signal = Data_total[N_start:N_end,2]- Data_total[nn_start, 2]*ones(length(t)); # sin (2π f t)

# Fourier Transform of it
F = fft(signal)/length(t)|>fftshift ;
time_domain = plot(t, signal, title = "Signal");
freqs = fftfreq(length(t), 1.0/Ts) |>fftshift;
freq_domain = plot(freqs, abs.(F), title = "Spectrum",xlim=(0, +30))
plot(time_domain, freq_domain, layout = (2,1))
	
end

# ╔═╡ 17aba62a-3ab3-11eb-11e6-2f95e58f627d
size(Data_total)

# ╔═╡ c8a11654-3ac7-11eb-344f-d55a714ae798
ind_FFT = findall(x-> x>0.0001, abs.(F))

# ╔═╡ 07aebc52-3ac8-11eb-2beb-17441fb89a9b
hcat(freqs[ind_FFT], abs.(F[ind_FFT]))

# ╔═╡ Cell order:
# ╠═0b1acabe-3ab6-11eb-3141-371834a04fbd
# ╠═97a4c138-3ab6-11eb-0365-d7c2e6a7ffe5
# ╠═50e7b37a-3ab6-11eb-184b-71e50ff72c3d
# ╠═afb55dae-3ab7-11eb-04fb-31cfab44f0bc
# ╠═adac1c30-3ac3-11eb-2d7a-c34e1278251f
# ╠═5cf4d328-3ab8-11eb-1948-5faeb828eca8
# ╠═071e8caa-3ab8-11eb-1e4f-b159d87cabe1
# ╠═50830bfc-3ab6-11eb-0c58-bfa15d7747d0
# ╠═5021084e-3ab6-11eb-15d3-279d2c4b11f0
# ╠═32cd8376-3ac0-11eb-2061-0f28dc1a9e77
# ╠═c8666f5a-3ac1-11eb-0025-95a05abc171d
# ╠═57c399e4-3ac1-11eb-3249-91634004453c
# ╠═7f4d6a52-3aa9-11eb-06c9-c14eb7146d4b
# ╠═df33fbb2-3ac1-11eb-1fe7-2f289db520b5
# ╠═3e88fd4c-3ac2-11eb-354d-f19f982f34aa
# ╠═49a76a08-3ac2-11eb-0043-d7fba49a40ea
# ╠═1a590fa6-3ac3-11eb-388f-1d68e504ad14
# ╠═81c43378-3ac3-11eb-02c1-8b3f25d8b640
# ╠═505751d2-3aab-11eb-1642-c55919d3ea8f
# ╠═6aa2a052-3aab-11eb-2420-7bb2200dfb25
# ╠═6c54fca0-3ac5-11eb-3c83-31e7a826392f
# ╠═c0d9925c-3ac5-11eb-2a23-fd1ab847886b
# ╠═bab146bc-3ac5-11eb-2b8d-ff2fe74d6ba3
# ╠═7502cebe-3aab-11eb-1318-e5aea2a6c9da
# ╠═37dce552-3ac6-11eb-19e0-cd99abf02a5c
# ╠═85babbf4-3aab-11eb-1418-2140dcd88435
# ╠═9b61420a-3ac6-11eb-0fd4-6f845ad884ae
# ╠═ac95f44e-3ac6-11eb-3a09-b77cb8f5e850
# ╠═4727a398-3ab0-11eb-35b6-5d4d0b4fb102
# ╠═8e3d5994-3aab-11eb-1daa-77db968d33c4
# ╠═17aba62a-3ab3-11eb-11e6-2f95e58f627d
# ╠═fe1ca76a-3aa8-11eb-1259-a1d5fdb20f4a
# ╠═c8a11654-3ac7-11eb-344f-d55a714ae798
# ╠═07aebc52-3ac8-11eb-2beb-17441fb89a9b
