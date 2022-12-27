import os
import math
import random

import matplotlib.pyplot as plt
from matplotlib import rcParams
from codeflaws_version_control import Qr_excel
import pandas as pd
import numpy as np
import matplotlib

version = {
    1: [
        (r'../report/codeflawf/sbfl+mbfl-1/sbfl.xlsx', 'sbfl'),
        (r'../report/codeflawf/sbfl+mbfl-1/mbfl.xlsx', 'mbfl'),
        (r'../report/codeflawf/sbfl+mbfl-1/random_hmbfl.xlsx', 'random \n hmbfl'),
        (r'../report/codeflawf/sbfl+mbfl-1/h-h.xlsx', 'high-high \n combin'),
        (r'../report/codeflawf/sbfl+mbfl-1/h-l.xlsx', 'high-low \n combin'),
        (r'../report/codeflawf/sbfl+mbfl-1/l-l.xlsx', 'low-low \n combin'),
    ],
    2: [
        (r'../report/codeflawf/sbfl+mbfl-5/sbfl.xlsx', 'sbfl'),
        (r'../report/codeflawf/sbfl+mbfl-5/mbfl.xlsx', 'mbfl'),
        (r'../report/codeflawf/sbfl+mbfl-2/hmbfl.xlsx', 'random \n hmbfl'),
        (r'../report/codeflawf/sbfl+mbfl-2/h-h.xlsx', 'high-high \n combin'),
        (r'../report/codeflawf/sbfl+mbfl-2/h-l.xlsx', 'high-low \n combin'),
        (r'../report/codeflawf/sbfl+mbfl-2/l-l.xlsx', 'low-low \n combin'),
    ],
    5: [
        (r'../report/codeflawf/sbfl+mbfl-5/sbfl.xlsx', 'sbfl'),
        (r'../report/codeflawf/sbfl+mbfl-5/mbfl.xlsx', 'mbfl'),
        (r'../report/codeflawf/sbfl+mbfl-5/hmbfl.xlsx', 'random \n hmbfl'),
        (r'../report/codeflawf/sbfl+mbfl-5/h-h.xlsx', 'high-high \n combin'),
        (r'../report/codeflawf/sbfl+mbfl-5/h-l.xlsx', 'high-low \n combin'),
        (r'../report/codeflawf/sbfl+mbfl-5/l-l.xlsx', 'low-low \n combin'),
    ]
}

excel = version[1]
exam_type = {
    'best': 8,
    'average': 9,
    'worst': 10,
}

sheet = 'muti'

font = {'family': 'Times New Roman'}

class Line:
    def __init__(self, x, y, path):
        plt.plot(x, y)
        plt.ylabel("Exam")
        plt.xlabel("HomNum/FomNum")
        # plt.title('')
        # plt.vlines(x, 0.294, y, linestyle="dashed")
        plt.xticks(x, x)
        # plt.xticks(x, x)
        plt.savefig(path)
        plt.show()
        plt.close()


class MutiBline:
    def __init__(self, row, col):
        self.f, self.ax = plt.subplots(row, col, figsize=(12, 8))
        self.maxrow = row
        self.maxcol = col

        plt.setp(self.ax.flat,
                 xlabel="% of statements that need to be examined",
                 ylabel="% of faulty version",
                 # fontdict=font,
                 )

    def draw(self, data, name, row, col):
        position = self.ax[row][col]
        # position.xlabel("% of statements that need to be examined", fontdict=font)
        # position.ylabel("% of faulty version", fontdict=font)

        color_total = ['r', 'orange', 'y', 'lawngreen',
                       'aqua', 'dodgerblue', 'violet', 'crimson', 'pink',
                       'gray', 'maroon', 'darkgreen']
        colorlist = color_total[:len(data.keys())-1] + ['black']
        for key, value in Bline().getBlineData(data).items():
            position.plot(value['x'], value['y'], color=colorlist.pop(0), label=key)
        position.legend(loc="best", prop=font)
        # position.xticks(position.xticks()[0][1:-1], fontproperties='Times New Roman')
        # position.yticks(position.yticks()[0][1:-1], fontproperties='Times New Roman')
        title = 'Fig.%s %s' % (str(self.maxcol*row + col+1), name)
        position.title = title

    def show(self, path, name):
        self.f.tight_layout()
        self.f.subplots_adjust(left=0.15, top=0.95)
        # plt.title = name
        # plt.savefig(path)
        plt.show()
        plt.close()



class Bline:
    def __init__(self):
        self.date = dict()
        self.title = ''
        self.xlabel = 'method'
        self.ylabel = 'exam'
        self.name = ''
        self.point = 0

    def getBlineData(self, date):
        blineDate = dict()

        for key, value in date.items():
            data = sorted(value, key=lambda x: x, reverse=True)
            x, y = [], []
            exam = data[0]
            while exam > data[-1]:
                x.append(exam*100)
                y.append((len(data) - data.index(exam))/len(data)*100)
                for i in range(data.index(exam), len(data)):
                    if data[i] < exam:
                        exam = data[i]
                        break
            x.append(exam*100)
            y.append((len(data) - data.index(exam))/len(data)*100)
            blineDate[key] = {
                'x': x,
                'y': y,
            }
        return blineDate

    def exam(self, path):
        color_total = ['r', 'orange', 'green', 'dodgerblue', 'y',
                     'aqua', 'pink', 'violet', 'crimson', 'lawngreen',
                     'gray', 'maroon', 'darkgreen']


        marker_total = [
            r'$\bigodot$', r'$\coprod$', r'$\Delta$', r'$\Vert$', 'D',
        ]
        linestyle_total = [
            ':', '--', '--', '-.', '-'
        ]
        colorlist = color_total[:len(self.date.keys())-1] + ['black']
        # plt.xlabel("line %")
        plt.rcParams.update({"font.size":15})
        plt.xlabel("% of statements that need to be examined", fontdict=font)
        plt.ylabel("% of faulty version", fontdict=font)
        for key, value in self.date.items():
            data = sorted(value, key=lambda x: x, reverse=True)
            x, y = [], []
            exam = data[0]
            while exam > data[-1]:
                x.append(exam*100)
                y.append((len(data) - data.index(exam))/len(data)*100)
                for i in range(data.index(exam), len(data)):
                    if data[i] < exam:
                        exam = data[i]
                        break
            x.append(exam*100)
            y.append((len(data) - data.index(exam))/len(data)*100)
            # print(key)
            line = plt.plot(x, y, color=colorlist.pop(0), label=key, marker=marker_total.pop(0), linestyle=linestyle_total.pop(0), markevery=10)
            # line = plt.plot(x, y, label=key, linestyle=linestyle_total.pop(0), marker=marker_total.pop(0), markevery=10)
            # plt.plot(x, y, color=colorlist.pop(0), label=key, linestyle=linestyle_total.pop(0), marker=marker_total.pop(0))
            # if key == 'Metallaxis':
            #     line.set_linewidth(2.0)
            if self.point > 0:
                for i_x in range(len(x)-1):
                    if min(x[i_x], x[i_x+1]) <= self.point*100 and self.point*100 <= max(x[i_x], x[i_x+1]):
                        maxx, minx = max(x[i_x], x[i_x+1]), min(x[i_x], x[i_x+1])
                        maxy, miny = max(y[i_x], y[i_x+1]), min(y[i_x], y[i_x+1])
                        iy = miny + (maxy-miny)/(maxx-minx)*(self.point*100-minx)
                        print('%s: %s-%s, %s-%s \n %s-%s' % (key, x[i_x], y[i_x], x[i_x+1], y[i_x+1], self.point*100, iy))
        plt.legend(loc="lower right", prop=font)
        if self.name == '':
            filepath = os.path.join(path, 'Bline-%s.pdf' % self.title)
        else:
            filepath = os.path.join(path, '%s.pdf' % self.name)
        a = plt.xticks()
        # plt.xticks(plt.xticks()[0][1:-1], fontproperties='Times New Roman')
        plt.xticks([i for i in range(0, 101, 10)], fontproperties='Times New Roman')
        # plt.yticks(plt.yticks()[0][1:-1], fontproperties='Times New Roman')
        plt.yticks([i for i in range(0, 101, 10)], fontproperties='Times New Roman')

        plt.hlines([i for i in range(0, 101, 10)], 0, 100, colors='gray', linestyles='dashed')


        plt.tight_layout()

        # plt.savefig(filepath, bbox_inches='tight', pad_inches=0)
        plt.savefig(filepath)
        plt.subplots_adjust(left=0.115, right=0.999, top=0.999, bottom=0.115)
        plt.show()
        plt.close()
        print(filepath)


class Box:
    def __init__(self):
        self.date = dict()

        self.title = ''
        self.xlabel = 'method'
        self.ylabel = 'exam'
        self.detail = False
        self.name = ''

    def exam(self, path):
        # 设置中文
        # plt.rcParams['font.sans-serif']=['SimHei','Times New Roman']
        plotdate = []
        plotlabel = []
        for key, value in self.date.items():
            plotdate.append(value)
            plotlabel.append(key)
        plt.boxplot(plotdate, labels=plotlabel)
        plt.ylim(0.0, 1.0)

        # data1 = [1,2,3,4]
        # data2 = [5,6,7,8]
        # plt.boxplot([data1, data2],labels=['male', 'female'])

        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        # plt.title('1')

        if self.detail:
            for i,x in enumerate(plotdate):
                sorteddata = sorted(x)
                Q1 = sorteddata[int(len(sorteddata)*0.75)]
                Q2 = sorteddata[int(len(sorteddata)*0.5)]
                Q3 = sorteddata[int(len(sorteddata)*0.25)]
                Qmax = Q1 + (Q1-Q3)*1.5
                Qmin = Q3 - (Q1-Q3)*1.5
                ytext = [Q3, Q2, Q1]
                ytick = -1
                for j in ytext:
                    if ytick == -1:
                        ytick = round(j, 2)
                    else:
                        if round(j, 2) - ytick < 0.1:
                            ytick += 0.1
                        else:
                            ytick = round(j, 2)
                    plt.text(i+1.3, ytick, round(j, 2))

        if self.name == '':
            filepath = os.path.join(path, 'Box-%s.png' % self.title)
        else:
            filepath = os.path.join(path, '%s.png' % self.name)
        plt.savefig(filepath)
        plt.show()
        print(filepath)


class Violin:
    def __init__(self):

        self.date = dict()
        self.title = ''
        self.xlabel = 'method'
        self.ylabel = 'exam'
        self.detail = False

    def exam(self, path):
        fig, axes = plt.subplots(figsize=(12, 5))
        axes.set_title(self.title)
        axes.set_xlabel(self.xlabel)
        axes.set_ylabel(self.ylabel)
        axes.yaxis.grid(True)
        plt.ylim(0.0, 1.0)

        all_data = []
        all_label = []
        for key, value in self.date.items():
            all_data.append(value)
            all_label.append(key)

        axes.violinplot(all_data, showmeans=False, showmedians=True)
        axes.set_xticks([y + 1 for y in range(len(all_data))], )

        plt.setp(axes, xticks=[y + 1 for y in range(len(all_data))], xticklabels=all_label,)

        filepath = os.path.join(path, 'violine-%s.png' % self.title)
        plt.savefig(filepath)
        plt.show()
        print(filepath)

class Hist:
    def __init__(self, data):
        matplotlib.rcParams['font.sans-serif']=['SimHei']   # 用黑体显示中文
        matplotlib.rcParams['axes.unicode_minus']=False     # 正常显示负号
        # 随机生成（10000,）服从正态分布的数据
        """
        绘制直方图
        data:必选参数，绘图数据
        bins:直方图的长条形数目，可选项，默认为10
        normed:是否将得到的直方图向量归一化，可选项，默认为0，代表不归一化，显示频数。normed=1，表示归一化，显示频率。
        facecolor:长条形的颜色
        edgecolor:长条形边框的颜色
        alpha:透明度
        """
        plt.hist(data, bins=40, normed=0, facecolor="blue", edgecolor="black", alpha=0.7)
        # 显示横轴标签
        plt.xlabel("区间")
        # 显示纵轴标签
        plt.ylabel("频数/频率")
        # 显示图标题
        plt.title("频数/频率分布直方图")
        plt.show()

def timespend():
    # 设置中文
    plt.rcParams['font.sans-serif']=['SimHei']

    plt.figure(figsize=(8, 6), dpi=80)
    # 再创建一个规格为 1 x 1 的子图
    plt.subplot(1, 1, 1)
    # 柱子总数
    # 包含每个柱子对应值的序列
    datay = []
    datax = []

    maxy = 0


    for j, strategy in enumerate(excel):
        datax.append(strategy[1])

        total_time = 0
        total_num = 0
        date_excel = Qr_excel().read(strategy[0], sheet)
        for i in range(len(date_excel)):
            if pd.isna(date_excel[i, 4]):
                continue
            if date_excel[i, 4] == 0:
                continue
            total_num += 1
            total_time += date_excel[i, 4]

        if j == 0:
            avetime = total_time/total_num/1000000
            plt.text(j, avetime, f"{avetime:.2e}")
        else:
            avetime = total_time/total_num
            plt.text(j, avetime, round(avetime, 2))
        datay.append(avetime)
        if avetime > maxy:
            maxy = avetime




    # 柱子的宽度
    width = 0.35
    # 绘制柱状图, 每根柱子的颜色为紫罗兰色
    p2 = plt.bar(np.arange(len(tuple(datay))), tuple(datay), width)
    # 添加纵横轴的刻度
    plt.xticks(np.arange(len(tuple(datay))), tuple(datax))
    plt.yticks(np.arange(0, 101, 10))
    plt.yticks(np.arange(0, (maxy//10+1)*10+1, (maxy//10+1)))


    # 设置横轴标签
    # 设置纵轴标签
    plt.ylabel('time spend(s)')
    plt.xlabel('methods')
    # 添加标题
    plt.title('average timespend')
    # 添加图例
    plt.legend(loc="upper right")
    plt.savefig('timespend.png', dpi=400)
    plt.show()


def draw():
    plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
    plt.title('显示中文标题')
    plt.xlabel("横坐标")
    plt.ylabel("纵坐标")

    marker_total = [
        r'$\bigodot$', r'$\coprod$', r'$\Delta$', r'$\Vert$', 'D',
    ]
    linestyle_total = [
        ':', '--', '--', '-.', '-'
    ]
    color_total = ['r', 'orange', 'green', 'pink', 'y',
                   'aqua', 'dodgerblue', 'violet', 'crimson', 'lawngreen',
                   'gray', 'maroon', 'darkgreen']
    # for i, mark in enumerate(marker_total):
    for i, color in enumerate(color_total):
    # for i, linestyle in enumerate(linestyle_total):
        list1=[x for x in range(10)]
        list2=[i for x in range(10)]
        plt.plot(list1, list2, label=color, color=color)
        # plt.plot(list1, list2, color='black', label=mark, marker=marker_total[i], linestyle=linestyle_total[i])
        # line = plt.plot(list1, list2, label=i, linestyle=linestyle, marker=',', markevery=10)
        # if i == 2:
        #     line[0].set_linewidth(2.0)
    plt.legend()
    # plt.grid()#添加网格
    plt.show()

if __name__ == '__main__':
    # timespend()
    # top()
    # exam()
    draw()
    # print(math.log(2.5+1))
    # Violin()