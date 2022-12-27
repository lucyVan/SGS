from codeflaws_version_control import Qr_excel
import pandas as pd
import numpy as np
import os,sys
import datetime,time
import mbfl.command as mbfl
import mbfl.mbfl_for as mbfl_for
import util, json

# 根据sbfl结果和随机hmbfl结果进行组合策略
def combination_strategy():
    sheet = 'muti'

    path_mbfl = r'./report/codeflawf/random_hmbfl.xlsx'
    date_mbfl = Qr_excel().read(path_mbfl, sheet)

    path_sbfl = r'../report/codeflawf/sbfl.xlsx'
    date_sbfl = Qr_excel().read(path_sbfl, sheet)

    path_hh = r'../report/c/h-h.xlsx'
    date_hh = Qr_excel().read(path_hh, sheet)
    path_hl = r'../report/c/h-l.xlsx'
    date_hl = Qr_excel().read(path_hl, sheet)
    path_ll = r'../report/c/l-l.xlsx'
    date_ll = Qr_excel().read(path_ll, sheet)
    strategy = {
        0: date_hh,
        1: date_hl,
        2: date_ll,
    }


    passnum = 0
    for i in range(len(date_sbfl)):
        print(i)
        if np.isnan(date_sbfl[i, 4]):
            continue
        out_list = {0: [], 1: [], 2: []}
        fom_list = {0: [], 1: [], 2: []}
        fom_time = {0: [], 1: [], 2: []}
        out_dic = {0: {}, 1: {}, 2: {}}

        try:
            Fault_Record = eval(date_sbfl[i, 1])
            linenum = date_sbfl[i, 3]
            sus = eval(date_sbfl[i, 5])


            hom_out = eval(date_mbfl[i, 11])
            homs = eval(date_mbfl[i, 12])
            or_list = eval(date_mbfl[i, 13])
            hom_time = eval(date_mbfl[i, 14])
            continue
        except:
            passnum += 1
            # print('读取错误')
            continue

        prehalf = sus[:len(sus)//2]
        dehalf = sus[len(sus)//2:]

# -----------------------------------------------------------------
        print('获取hom')
        for j, hom_str in enumerate(homs):
            hom = eval(hom_str)
            strategy_type = 0
            # 0: h-h 1:h-l  2:l-l
            for fom in hom:
                if not any(fom[0] in sus_single for sus_single in prehalf):
                    strategy_type += 1
            out_list[strategy_type].append(hom_out[j])
            fom_list[strategy_type].append(hom)
            fom_time[strategy_type].append(hom_time[j])

            for fom in hom:
                line = fom[0]
                if line not in out_dic[strategy_type]:
                    out_dic[strategy_type][line] = []
                out_dic[strategy_type][line].append(hom_out[j])


# -----------------------------------------------------------------
        print('生成怀疑度')
        for j in range(3):
            touple = mbfl.get_toule(or_list, out_dic[strategy_type])
            sus_list = []
            for line in touple['f&p']:
                touple_list = touple['tf'], touple['f&p'][line]['f'], touple['f&p'][line]['p'], touple['f2p'], touple[
                    'p2f'],
                # sus = mbfl_for.metallaxis(touple_list)
                sus = mbfl_for.muse(touple_list)
                sus_list.append([line, sus])
            sus_list_sort = sorted(sus_list, key=lambda x: x[1], reverse=True)
            ranks = util.Evaluation().rank(sus_list_sort, Fault_Record)
            best, average, worst = [], [], []
            for fault in ranks:
                best.append(fault[1])
                average.append(fault[2])
                worst.append(fault[3])

# -----------------------------------------------------------------
            strategy.get(j)[i, 3] = linenum
            strategy.get(j)[i, 5] = sus_list_sort
            strategy.get(j)[i, 8] = best
            strategy.get(j)[i, 9] = average
            strategy.get(j)[i, 10] = worst
            strategy.get(j)[i, 4] = sum(fom_time[j])/1000000
            strategy.get(j)[i, 11] = out_list[j]
            strategy.get(j)[i, 12] = fom_list[j]
            strategy.get(j)[i, 13] = or_list
            strategy.get(j)[i, 14] = fom_time[j]

# -----------------------------------------------------------------
        print('存储数据')
        # date_mbfl[i, 4] = sum(eval(hom_time))
        Qr_excel().save(path_hh, sheet, strategy.get(0))
        Qr_excel().save(path_hl, sheet, strategy.get(1))
        Qr_excel().save(path_ll, sheet, strategy.get(2))
        # Qr_excel().save(path_mbfl, sheet, date_mbfl)
    print(passnum)


def combination_json_strategy():
    sheet = 'muti'

    filename = '../report/codeflawf/sbfl+mbfl-2/random-hmbfl-2.json'
    with open(filename) as f_obj:
        date_mbfl = json.load(f_obj)

    path_sbfl = r'../report/codeflawf/sbfl+mbfl-2/sbfl.xlsx'
    date_sbfl = Qr_excel().read(path_sbfl, sheet)


    path_hh = r'../report/codeflawf/sbfl+mbfl-2/h-h.xlsx'
    date_hh = Qr_excel().read(path_hh, sheet)
    path_hl = r'../report/codeflawf/sbfl+mbfl-2/h-l.xlsx'
    date_hl = Qr_excel().read(path_hl, sheet)
    path_ll = r'../report/codeflawf/sbfl+mbfl-2/l-l.xlsx'
    date_ll = Qr_excel().read(path_ll, sheet)
    path_hmbfl = r'../report/codeflawf/sbfl+mbfl-2/hmbfl.xlsx'
    date_hmbfl = Qr_excel().read(path_hmbfl, sheet)

    strategy = {
            0: date_hh,
            1: date_hl,
            2: date_ll,
            3: date_hmbfl,
        }

    for i in range(len(date_sbfl)):
        doc = date_sbfl[i, 0]
        if doc not in date_mbfl:
            continue
        if np.isnan(date_sbfl[i, 4]):
            continue

        try:
            Fault_Record = eval(date_sbfl[i, 1])    # 故障行
            linenum = date_sbfl[i, 3]    # 源程序代码行数
            sus = eval(date_sbfl[i, 5])    # sbfl下怀疑度列表

            or_list = date_mbfl[doc]['sorce']
            fomnum = date_mbfl[doc]['fomnum']    # 一阶变异体数量
            homs = date_mbfl[doc]['homs']    # 二阶变异体信息

        except Exception as e:
            continue
        print(i, doc)


        out_list = {0: [], 1: [], 2: [], 3: []}
        fom_time = {0: [], 1: [], 2: [], 3: []}
        fom_list = {0: [], 1: [], 2: [], 3: []}
        hom_list = {0: [], 1: [], 2: [], 3: []}
        out_dic = {0: {}, 1: {}, 2: {}, 3: {}}
        rank_best = {0: [], 1: [], 2: [], 3: []}
        rank_average = {0: [], 1: [], 2: [], 3: []}
        rank_worst = {0: [], 1: [], 2: [], 3: []}

        prehalf = sus[:len(sus)//2]
        dehalf = sus[len(sus)//2:]

        # -----------------------------------------------------------------
        print('获取hom')
        for j, hom in enumerate(homs):
            strategy_type = 0
            # 0: h-h 1:h-l  2:l-l
            out = hom['out']
            time = hom['time']
            message = eval(hom['message'])


            # 跳过在同一行生成fom的hom
            breakthishom = False
            for fom1 in message:
                for fom2 in message:
                    if fom1[0] == fom2[0] and fom1[2] != fom2[2]:
                        breakthishom = True
            if breakthishom:
                continue

            for fom in message:
                if not any(fom[0] in sus_single for sus_single in prehalf):
                    strategy_type += 1
            out_list[strategy_type].append(out)
            hom_list[strategy_type].append(message)
            fom_time[strategy_type].append(time)

            out_list[3].append(out)
            hom_list[3].append(message)
            fom_time[3].append(time)

            for fom in message:
                line = fom[0]
                if line not in out_dic[strategy_type]:
                    out_dic[strategy_type][line] = []
                if line not in out_dic[3]:
                    out_dic[3][line] = []
                out_dic[strategy_type][line].append(out)
                out_dic[3][line].append(out)


        # -----------------------------------------------------------------
        print('生成怀疑度')
        for j in range(4):
            touple = mbfl.get_toule(or_list, out_dic[j])
            sus_list = []
            for line in touple['f&p']:
                touple_list = touple['tf'], touple['f&p'][line]['f'], touple['f&p'][line]['p'], touple['f2p'], touple[
                    'p2f'],
                # sus = mbfl_for.metallaxis(touple_list)
                sus = mbfl_for.muse(touple_list)
                sus_list.append([line, sus])
            sus_list_sort = sorted(sus_list, key=lambda x: x[1], reverse=True)
            ranks = util.Evaluation().rank(sus_list_sort, Fault_Record)
            best, average, worst = [], [], []
            for fault in ranks:
                best.append(fault[1])
                average.append(fault[2])
                worst.append(fault[3])

            # -----------------------------------------------------------------
            strategy.get(j)[i, 3] = linenum
            strategy.get(j)[i, 5] = sus_list_sort
            strategy.get(j)[i, 8] = best
            strategy.get(j)[i, 9] = average
            strategy.get(j)[i, 10] = worst
            strategy.get(j)[i, 4] = sum(fom_time[j])/1000000
            strategy.get(j)[i, 11] = out_list[j]
            strategy.get(j)[i, 12] = fom_list[j]
            strategy.get(j)[i, 13] = or_list
            strategy.get(j)[i, 14] = fom_time[j]

        # -----------------------------------------------------------------
        print('存储数据')
        # date_mbfl[i, 4] = sum(eval(hom_time))
        Qr_excel().save(path_hh, sheet, strategy.get(0))
        Qr_excel().save(path_hl, sheet, strategy.get(1))
        Qr_excel().save(path_ll, sheet, strategy.get(2))
        Qr_excel().save(path_hmbfl, sheet, strategy.get(3))




if __name__ == '__main__':
    combination_json_strategy()
    # t1 = datetime.datetime.now()
    # time.sleep(0.3)
    # t2 = datetime.datetime.now()
    # t3 = (t2-t1).microseconds
    # pass


