#f = open("./path_sdrive.txt")
f = open("./freq41Hz/t3s2_type4_inj.inf"); # in MacOS, path use ./. IN windows, path use \\

lines = readlines(f)

# find the column names
names = Array{String}(undef,length(lines));

for i=1:length(lines)
	ind_desc = findfirst("Desc=\"", lines[i]); 
        name_str = ind_desc[end]+1;
        ind_desc1 = findnext("\"", lines[i], name_str);
        name_end = ind_desc1[1]-1;
        names[i] = lines[i][name_str:name_end];
end

# compute how many files will be used to store the data
N_files = Int(floor(length(names)/10)+1);

filenames = ["./freq41Hz/t3s2_type4_inj_"*string(i, pad=2)*".out" for i=1:N_files]

names

#using CSV
#CSV.File(filenames[1])

using DelimitedFiles

Data_total,heading =readdlm(filenames[1],  header=true)
size(Data_total)

#for i=1:N_files
#open(filenames[1]) do io
#open(filenames[1]*"_w","w") do io1
#        k =1; 
#  while true
#         
#       line =  readline(io)
#       line1 = line[4:end]
#            if(k>1)
#               write(io1, line1)
#            end
#            k = k+1;
#  end
#end
#end

        
#end

size(Data_total)

for i =2: N_files
    data, heading =readdlm(filenames[i],  header=true)
    Data_total = [Data_total data[:,2:end]]
end

size(Data_total)

using Plots
Using GKS
plot(Data_total[:,1], Data_total[:,2:20], title ="41 Hz perturbation")

# time interval is 0.0005 seconds or 2000 Hz

using FFTW
using AbstractFFTs


tt = Data_total[:,1];
Ts = tt[2]-tt[1];
# Start time
t0 = 0.9; 
N = 8000-1;
N_start = Int(floor(t0/Ts));
tmax = t0 + N * Ts
N_end = N_start + N;


# time coordinate
#t = [t0:Ts:tmax;]
t = Data_total[N_start:N_end,1]

# signal
signal = Data_total[N_start:N_end,2]  # sin (2π f t)

# Fourier Transform of it
F = fft(signal)
time_domain = plot(t, signal, title = "Signal")
freqs = fftfreq(length(t), 1.0/Ts)
freq_domain = plot(freqs, abs.(F), title = "Spectrum",xlim=(-100, +100))
plot(time_domain, freq_domain, layout = 2)

abs.(F)


size(t)






