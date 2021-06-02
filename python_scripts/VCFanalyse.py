from pandas import DataFrame


def read_vcf(fh):
	metadata = []
	DF = None

	for i in fh:
		if(i[0:2]=='##'):
			metadata.append(i.strip())
		elif(i[0:6]=='#CHROM'):
			DF = DataFrame( columns=i.strip().split() )
		else:
			#DF.append( { DF.columns[numj]:j for numj,j in enumerate(i.split()) } , ignore_index=True)
			try:
				DF = DF.append( DataFrame( [[j for j in i.strip().split()]] , columns= list(DF.columns) ) )
			except:
				new_columns = [ 'col'+str(j+1) for j in range(len(DF.columns), len(DF.columns)+len(i.split())-len(DF.columns) ) ]
				for j in new_columns: DF[j] = None
				DF = DF.append( DataFrame( [[j for j in i.strip().split()]] , columns= list(DF.columns) ) )

	return DF, metadata


def main():
	#example
	#fh = open('NA19461.final.manta.svtyped.vcf')
	fh = open('NA19461.final.breakdancer.vcf')
	records , metadata = read_vcf(fh)

	print(records)

	


if(__name__=='__main__'):
	main()