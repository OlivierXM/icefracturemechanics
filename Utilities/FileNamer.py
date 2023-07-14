import datetime
import os

def FileNamer(arg1,arg2 = "Test"):
    if (not isinstance(arg1,str)):
        raise Exception('Arg1', 'must be str')

    curDate = datetime.datetime.now()
    curYear = curDate.date().year - 2000 # Sorry 20th century
    curMonth = curDate.date().month
    curDay = curDate.date().day
    baseStr = '{:02.0f}{:02.0f}{:02.0f}'.format(curYear,curMonth,curDay)
    numDate = curYear*10000 + curMonth*100 + curDay # YYMMDD (i.e. 210607)

    if (not os.path.isfile("../SaveLog.txt")):
        print("SaveLog.txt not found, creating...")
        fNew = open("../SaveLog.txt","w")
        fNew.write("{0}\n".format(str(numDate)))
        fNew.write("00\n")
        fNew.close()

    f = open('../SaveLog.txt','r')
    testDate = f.readline()
    testNumber = f.readline()
    if(int(testDate[:-1]) < numDate):
        f.close()
        f = open('../SaveLog.txt','w')
        f.write(str(numDate)+"\n")
        f.write("01\n")
        fileNumber = "00\n"
        f.close()
    else:
        fileNumber = testNumber
        f.close()
        f = open('../SaveLog.txt','w')
        f.write(testDate)
        f.write("{:02.0f}\n".format(int(testNumber)+1))
        f.close()

    fileName = baseStr + "_" + arg2 + "_" + fileNumber[:-1] + "_" + arg1

    return fileName

def FileNamerV2(arg1,arg2 = "Test", num=0):
    curDate = datetime.datetime.now()
    curYear = curDate.date().year - 2000 # Sorry 20th century
    curMonth = curDate.date().month
    curDay = curDate.date().day
    baseStr = '{:02.0f}{:02.0f}{:02.0f}'.format(curYear,curMonth,curDay)
    numDate = curYear*10000 + curMonth*100 + curDay # YYMMDD (i.e. 210607)

    fileNumber = "{:02.0f}".format(int(num))

    fileName = baseStr + "_" + arg2 + "_" + fileNumber + "_" + arg1

    return fileName

def ReturnNum():
    curDate = datetime.datetime.now()
    curYear = curDate.date().year - 2000 # Sorry 20th century
    curMonth = curDate.date().month
    curDay = curDate.date().day
    baseStr = '{:02.0f}{:02.0f}{:02.0f}'.format(curYear,curMonth,curDay)
    numDate = curYear*10000 + curMonth*100 + curDay # YYMMDD (i.e. 210607)

    if (not os.path.isfile("../SaveLog.txt")):
        print("SaveLog.txt not found, creating...")
        fNew = open("../SaveLog.txt","w")
        fNew.write("{0}\n".format(str(numDate)))
        fNew.write("00\n")
        fNew.close()

    f = open('../SaveLog.txt','r')
    testDate = f.readline()
    testNumber = f.readline()
    if(int(testDate[:-1]) < numDate):
        f.close()
        f = open('../SaveLog.txt','w')
        f.write(str(numDate)+"\n")
        f.write("01\n")
        fileNumber = "00\n"
        f.close()
    else:
        fileNumber = testNumber
        f.close()
        f = open('../SaveLog.txt','w')
        f.write(testDate)
        f.write("{:02.0f}\n".format(int(testNumber)+1))
        f.close()

    print(fileNumber)
    return 0
