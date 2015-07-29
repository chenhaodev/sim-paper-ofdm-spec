import os
import re
import sys
filename=sys.argv[1]
nfilename=sys.argv[1]+".tmp"
f=open(filename, "r")
lines=f.readlines()
fn=open(nfilename, "w")

pattern = 'IMAGE='
end_node = '/>'
#pre_text1='<richcontent TYPE="NODE"><html><head></head><body><img src='
pre_text1='<richcontent TYPE="NODE"><html><head/><body><img src='
pos_text1='/></body></html></richcontent></node>'

pre_text2='<richcontent TYPE="NODE"><html><head/><body><img src='
pos_text2='/></body></html></richcontent>'

for i in range(0,len(lines)):
	line=lines[i]
	if((line.find(pattern)!=-1)):
		text = line.split(pattern)[1].strip()
		text = text.replace(">", "")
		text = text.replace("<", "")
		src = re.split(r"\s",text)[0]
		if((line.find(end_node)!=-1)):
			line = line.replace("/>", ">")
			new_line = pre_text1 + src + pos_text1
		else:
			new_line = pre_text2 + src + pos_text2
		fn.write(line)
		fn.write(new_line + '\n')
	else:
		fn.write(line)
cmd1="mv " + nfilename + " " + filename
os.system(cmd1)
