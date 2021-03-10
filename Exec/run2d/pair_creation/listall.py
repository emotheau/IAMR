from glob import glob
import os.path

dirslist =[int(os.path.basename(x)[5:]) for x in glob("/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/data/256/*")]
#print(dirslist.shape)
print(max(dirslist))
print(len(dirslist))
