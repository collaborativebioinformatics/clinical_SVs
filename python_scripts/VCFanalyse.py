from pandas import DataFrame


def read_vcf(fh):
	metadata = []
	DF = None

	for i in fh:
		if(i[0:2]=='##'):
			metadata.append(i)
		elif(i[0:6]=='#CHROM'):
			DF = DataFrame( columns=i.split() )
		else:
			#DF.append( { DF.columns[numj]:j for numj,j in enumerate(i.split()) } , ignore_index=True)
			DF = DF.append( DataFrame( [[j for j in i.split()]] , columns= list(DF.columns) ) )
	
	return DF, metadata


def main():
	fh = open('NA19461.final.manta.svtyped.vcf')
	data = read_vcf(fh)


if(__name__=='__main__'):
	main()