import gbm
import os

yy = ['08','09','10','11','12','13','14','15','16','17','18','19']
mm = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
dd = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31']

for i in range(len(yy)):
    for j in range(len(mm)):
        if int(yy[i]) % 4 != 0 and int(mm[j]) == 2:  # Check for leap year
            for k in range(len(dd)-2):
                link = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/20"+yy[i]+"/"+mm[j]+"/"+dd[k]+"/current/glg_poshist_all_"+yy[i]+mm[j]+dd[k]+"_v00.fit"

                dest = "/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/glg_poshist_all_"+yy[i]+mm[j]+dd[k]+"_v00.fit"

                os.system("curl  %s > %s" % (link, dest))
        if int(yy[i]) == 4 or int(yy[i]) == 6 or int(yy[i]) == 9 or int(yy[i]) == 11:  # Check for months with 30 days
            for k in range(len(dd)-1):
                link = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/20"+yy[i]+"/"+mm[j]+"/"+dd[k]+"/current/glg_poshist_all_"+yy[i]+mm[j]+dd[k]+"_v00.fit"

                dest = "/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/glg_poshist_all_"+yy[i]+mm[j]+dd[k]+"_v00.fit"

                os.system("curl  %s > %s" % (link, dest))
        else:
            for k in range(len(dd)):
                link = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/20"+yy[i]+"/"+mm[j]+"/"+dd[k]+"/current/glg_poshist_all_"+yy[i]+mm[j]+dd[k]+"_v00.fit"

                dest = "/Users/RonnyNguyen/Desktop/FermiGBMSummer2020/data_dir/glg_poshist_all_"+yy[i]+mm[j]+dd[k]+"_v00.fit"

                os.system("curl  %s > %s" % (link, dest))
'''
"https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/YYYY/MM/DD/current/glg_poshist_all_YYMMDD_v00.fit"
'''
