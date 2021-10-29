import matplotlib.pyplot as plt

depth_file = "amostra-lbb_depth.sample_summary"
dic = {}
f = open(depth_file)
for line in f:
	break
for line in f:
	fields = line.split(',')
	sample = fields[0]
	d5 = float(fields[6])
	d10 = float(fields[7])
	d20 = float(fields[8])
	d50 = float(fields[9])
	d100 = float(fields[10])
	d150 = float(fields[11])
	d250 = float(fields[12])
	if sample not in dic:
		dic[sample] = [d5,d10,d20,d50,d100,d150,d250]
	break

y = [5,10,20,50,100,150,250]
for k,v in dic.items():
	x = v
	plt.plot(y,x,label=k)
plt.ylim([0,100])
plt.axvline(x=5, color='black', linestyle='--')
plt.axvline(x=10, color='black', linestyle='--')
plt.axvline(x=20, color='black', linestyle='--')
plt.axvline(x=50, color='black', linestyle='--')
plt.axvline(x=100, color='black', linestyle='--')
plt.axvline(x=150, color='black', linestyle='--')
plt.axvline(x=250, color='black', linestyle='--')

plt.axhline(y=50, color='r', linestyle='--')

plt.legend(loc='center left',ncol=5, bbox_to_anchor=(1, 0.5))
plt.xlabel('Depth of Coverage') 
plt.ylabel('% on-target bases covered') 
plt.xticks(y)
plt.show()

