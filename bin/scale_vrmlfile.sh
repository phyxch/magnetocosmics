echo "3a \\
Transform {  
3a \\
scale 0.0000001 0.0000001 0.0000001 
3a \\
children [ 
/End of file/ i \\
]} 
s/Shape/,Shape/g 
s/Viewpoint/,Viewpoint/g" >file_script1
sed --file=file_script1  $1 > $2
rm file_script1


