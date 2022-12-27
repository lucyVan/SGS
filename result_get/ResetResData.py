import os
import json
import datetime
import json
import shutil
# from codeflaws_version_control import Qr_excel


def GetKilList(orout, mutout):
    killlist = []
    for i in range(len(orout)):
        if orout[i] == mutout[i]:
            killlist.append(1)
        else:
            killlist.append(0)
    return killlist




def ResetSbflData():
    sbflpath = os.path.join(os.getcwd(), '../report/SBFL')
    chmbflpath = os.path.join(os.getcwd(), '../report/CHMBFL/mutinfo-single')
    for file in os.listdir(chmbflpath):

        chmbfl_file_path = os.path.join(chmbflpath, file)
        print('%s 读取 %s' % (datetime.datetime.now(), chmbfl_file_path), end='')
        try:
            with open(chmbfl_file_path) as f_obj:
                data_chmbfl = json.load(f_obj)
        except Exception as e:
            print('\r %s 读取失败 %s %s' % (datetime.datetime.now(), chmbfl_file_path, e))
            continue
        doc = list(data_chmbfl.keys())[0]

        sbfl_file_path = os.path.join(sbflpath, 'data_%s.json' % doc)
        try:
            with open(sbfl_file_path) as f_obj:
                data_sbfl = json.load(f_obj)
        except:
            print('\r %s 读取失败 %s' % (datetime.datetime.now(), sbfl_file_path))
            continue

        sbfldata = data_sbfl[doc]['sbfl']
        data_chmbfl[doc]['sbfl'] = sbfldata

        with open(chmbfl_file_path, 'w') as f_obj:
            json.dump(data_chmbfl, f_obj)
        print('\r %s 存储完成 %s' % (datetime.datetime.now(), chmbfl_file_path))
    return



def delDataLackVersion():
    readpath = os.path.join(os.getcwd(), '../report/CHMBFL/mutinfo')
    # readpath = os.path.join(os.getcwd(), '../report/CHMBFL/mutinfo-needadd')
    for i, file in enumerate(os.listdir(readpath)):
        filereadpath = os.path.join(readpath, file)

        print('%s 读取 %s - %s' % (datetime.datetime.now(), file, i), end='')
        try:
            with open(filereadpath) as f_obj:
                data_json = json.load(f_obj)
        except:
            print('\r %s %s - %s 读取失败' % (datetime.datetime.now(), file, i))
            continue

        doc = list(data_json.keys())[0]
        data = data_json[doc]
        del data_json

        try:
            if len(data['fom_list']) == 0:
                print('\r %s %s - %s 缺少fom数据' % (datetime.datetime.now(), file, i))
                continue
        except:
            print('\r %s %s - %s 缺少fom数据' % (datetime.datetime.now(), file, i))
            continue

        try:
            if data['homnum'] == 0:
                print('\r %s %s - %s 缺少hom数据' % (datetime.datetime.now(), file, i))
                continue
        except:
            print('\r %s %s - %s 缺少hom数据' % (datetime.datetime.now(), file, i))
            continue
        print('')

    return

def moveData():
    readpath = os.path.join(os.getcwd(), '../report/CHMBFL/1/NewData')
    targetpath = os.path.join(os.getcwd(), '../report/CHMBFL/NewData')
    for i, file in enumerate(os.listdir(readpath)):
        filereadpath = os.path.join(readpath, file)
        filetargetpath = os.path.join(targetpath, file)


        if not os.path.exists(filetargetpath):
            shutil.copy(filereadpath, filetargetpath)
            print('\r %s %s - %s 移动完成' % (datetime.datetime.now(), file, i))
        else:
            print('\r %s %s - %s 文件存在' % (datetime.datetime.now(), file, i))

    return


def A():
    sbfldatapath = '../report/SBFL'
    vc_path = '../report/data_codeflaws.xlsx'
    sheet = 'muti'
    date_read = Qr_excel().read(vc_path, sheet)
    for i in range(len(date_read)):
        Fault_Record = eval(date_read[i, 1])
        doc = date_read[i, 0]
        if not 'data_%s.json'%doc in os.listdir(sbfldatapath):
            continue
        sbflpath = os.path.join(sbfldatapath, 'data_%s.json' % doc)

        with open(sbflpath) as f_obj:
            data_json = json.load(f_obj)
        data_json[doc]['Fault_Record'] = Fault_Record

        with open(sbflpath, 'w') as f_obj:
            json.dump(data_json, f_obj)
        print('存储完成 %s' % sbflpath)











if __name__ == '__main__':
    # print(len(os.listdir('../report/CHMBFL/1/mutinfo')))
    # ResetCHMBFLData()
    ResetSbflData()
    # delDataLackVersion()
    # moveData()
    # A()
