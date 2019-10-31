import os

string_name1='s8a20130408_00016_00'    #==========================================HERE:a,b,c,d
string_name2='s8b20130408_00016_00'    #==========================================HERE:a,b,c,d
string_name3='s8c20130408_00016_00'    #==========================================HERE:a,b,c,d
string_name4='s8d20130408_00016_00'    #==========================================HERE:a,b,c,d

for i in range(1, 53):      #==========================================HERE:depending on the file numbers
    print i
    if i%2==0: #even
        try:
            os.rename('./'+string_name1+'{:02d}.sdf'.format(i), './even_folder/'+string_name1+'{:02d}.sdf'.format(i))
        except:
            pass
        try:
            os.rename('./'+string_name2+'{:02d}.sdf'.format(i), './even_folder/'+string_name2+'{:02d}.sdf'.format(i))
        except:
            pass
        try:
            os.rename('./'+string_name3+'{:02d}.sdf'.format(i), './even_folder/'+string_name3+'{:02d}.sdf'.format(i))
        except:
            pass
        try:
            os.rename('./'+string_name4+'{:02d}.sdf'.format(i), './even_folder/'+string_name4+'{:02d}.sdf'.format(i))
        except:
            pass
    else:
        try:
            os.rename('./'+string_name1+'{:02d}.sdf'.format(i), './odd_folder/'+string_name1+'{:02d}.sdf'.format(i))
        except:
            pass
        try:
            os.rename('./'+string_name2+'{:02d}.sdf'.format(i), './odd_folder/'+string_name2+'{:02d}.sdf'.format(i))
        except:
            pass
        try:
            os.rename('./'+string_name3+'{:02d}.sdf'.format(i), './odd_folder/'+string_name3+'{:02d}.sdf'.format(i))
        except:
            pass
        try:
            os.rename('./'+string_name4+'{:02d}.sdf'.format(i), './odd_folder/'+string_name4+'{:02d}.sdf'.format(i))
        except:
            pass

#{0:02d} for python 2.6
#{:02d} for python 2.7
