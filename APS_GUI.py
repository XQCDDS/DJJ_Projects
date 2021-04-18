import wx
import APS_3d
import APS_2d
import re

DNA_SIZE = 24  # 基因长度
POP_SIZE = 200  # 种群个数
CROSSOVER_RATE = 0.8  # 交叉概率
MUTATION_RATE = 0.005  # 变异概率
N_GENERATIONS = 50  # 迭代次数
X_BOUND = [-3, 3]  # x取值范围
Y_BOUND = [-3, 3]  # y取值范围
Z_BOUND = [-10, 10]  # y取值范围
ACCURACY = 100  # 取点精度
EXPRESSION_2D = ""  # 2D表达式
EXPRESSION_3D = ""  # 3D表达式
CURVE_POINT_2D = "(1,2),(2,1),(3,0.3)"  # 2D曲线拟合点
CURVE_POINT_3D = "(2,2,1),(2,0.3,0.15),(2,1,0.5),\
    (1,1,1),(1,2,2),(1,0.3,0.3),\
    (3,2,0.67),(3,1,0.33),(3,0.3,0.1)"  # 3D曲线拟合点


# 自定义窗口类MyFrame
class APSFrame(wx.Frame):
    def __init__(self, superior):
        wx.Frame.__init__(self, parent=superior, title=u'APS_GUI界面')
        # self.Centre()
        self.SetSize(800, 700)

        panel = wx.Panel(self, -1)
        # 添加静态文本控件
        label_dna_size = wx.StaticText(parent=panel, id=101, label=u"基因长度 ：")
        label_pop_size = wx.StaticText(parent=panel, id=102, label=u"种群个数 ：")
        label_crossover_rate = wx.StaticText(parent=panel, id=103, label=u"交叉概率 ：")
        label_mutation_rate = wx.StaticText(parent=panel, id=104, label=u'变异概率 ：')
        label_n_generations = wx.StaticText(parent=panel, id=105, label=u'世代次数 ：')
        label_x_bound = wx.StaticText(parent=panel, id=106, label=u'x取值范围 ：')
        label_y_bound = wx.StaticText(parent=panel, id=107, label=u'y取值范围 ：')
        label_z_bound = wx.StaticText(parent=panel, id=108, label=u'z取值范围 ：')
        label_accuracy = wx.StaticText(parent=panel, id=109, label=u'取点精度 ：')
        label_expression_2d = wx.StaticText(parent=panel, id=110, label=u'2D表达式 ：')
        label_expression_3d = wx.StaticText(parent=panel, id=111, label=u'3D表达式 ：')
        label_curve_point_2d = wx.StaticText(parent=panel, id=112, label=u'2D真实数据点 ：')
        label_curve_point_3d = wx.StaticText(parent=panel, id=113, label=u'3D真实数据点 ：')
        # 添加文本框
        self.input_dna_size = wx.TextCtrl(parent=panel, id=201, value=str(DNA_SIZE))
        self.input_pop_size = wx.TextCtrl(parent=panel, id=202, value=str(POP_SIZE))
        self.input_crossover_rate = wx.TextCtrl(parent=panel, id=203, value=str(CROSSOVER_RATE))
        self.input_mutation_rate = wx.TextCtrl(parent=panel, id=204, value=str(MUTATION_RATE))
        self.input_n_generations = wx.TextCtrl(parent=panel, id=205, value=str(N_GENERATIONS))
        self.input_x_bound = wx.TextCtrl(parent=panel, id=206, value=str(X_BOUND))
        self.input_y_bound = wx.TextCtrl(parent=panel, id=207, value=str(Y_BOUND))
        self.input_z_bound = wx.TextCtrl(parent=panel, id=208, value=str(Z_BOUND))
        self.input_accuracy = wx.TextCtrl(parent=panel, id=209, value=str(ACCURACY))
        self.input_expression_2d = wx.TextCtrl(parent=panel, id=210, value=str(EXPRESSION_2D)
                                               , size=(100, 100), style=wx.TE_MULTILINE)
        self.input_expression_3d = wx.TextCtrl(parent=panel, id=211, value=str(EXPRESSION_3D)
                                               , size=(100, 100), style=wx.TE_MULTILINE)
        self.input_curve_point_2d = wx.TextCtrl(parent=panel, id=212, value=str(CURVE_POINT_2D)
                                                , size=(100, 100), style=wx.TE_MULTILINE)
        self.input_curve_point_3d = wx.TextCtrl(parent=panel, id=213, value=str(CURVE_POINT_3D)
                                                , size=(100, 100), style=wx.TE_MULTILINE)

        # 添加按钮
        self.button_APS_3d = wx.Button(parent=panel, id=1, label=u'三维APS排程')
        self.button_APS_2d = wx.Button(parent=panel, id=2, label=u'二维APS排程')
        self.button_Clear = wx.Button(parent=panel, id=3, label=u'初始化')
        self.button_Back = wx.Button(parent=panel, id=4, label=u'返回')

        # 绑定事件处理函数
        self.Bind(wx.EVT_BUTTON, handler=self.OnButton_APS_3d, source=self.button_APS_3d)
        self.Bind(wx.EVT_BUTTON, handler=self.OnButton_APS_2d, source=self.button_APS_2d)
        self.Bind(wx.EVT_BUTTON, handler=self.OnButton_Clear, source=self.button_Clear)
        self.Bind(wx.EVT_BUTTON, self.OnButton_Back, source=self.button_Back)

        self.panel = panel

        # 创建水平方向box布局管理器（默认水平方向）
        hbox1 = wx.BoxSizer()
        # 将两个button添加到hbox布局管理器中
        hbox1.Add(label_dna_size, 101, wx.ALIGN_CENTER | wx.LEFT, border=10)
        hbox1.Add(self.input_dna_size, 201, wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=10)
        hbox1.Add(label_pop_size, 102, wx.ALIGN_CENTER | wx.LEFT, border=10)
        hbox1.Add(self.input_pop_size, 202, wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=10)
        hbox1.Add(label_accuracy, 109, wx.ALIGN_CENTER | wx.LEFT, border=10)
        hbox1.Add(self.input_accuracy, 209, wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=10)

        hbox2 = wx.BoxSizer()
        hbox2.Add(label_x_bound, 106, wx.ALIGN_CENTER | wx.LEFT, border=10)
        hbox2.Add(self.input_x_bound, 206, wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=10)
        hbox2.Add(label_y_bound, 107, wx.ALIGN_CENTER | wx.LEFT, border=10)
        hbox2.Add(self.input_y_bound, 207, wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=10)
        hbox2.Add(label_z_bound, 108, wx.ALIGN_CENTER | wx.LEFT, border=10)
        hbox2.Add(self.input_z_bound, 208, wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=10)

        hbox3 = wx.BoxSizer()
        hbox3.Add(label_crossover_rate, 103, wx.ALIGN_CENTER | wx.LEFT, border=10)
        hbox3.Add(self.input_crossover_rate, 203, wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=10)
        hbox3.Add(label_mutation_rate, 104, wx.ALIGN_CENTER | wx.LEFT, border=10)
        hbox3.Add(self.input_mutation_rate, 204, wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=10)
        hbox3.Add(label_n_generations, 105, wx.ALIGN_CENTER | wx.LEFT, border=10)
        hbox3.Add(self.input_n_generations, 205, wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=10)

        hbox4 = wx.BoxSizer()
        hbox4.Add(label_expression_3d, 111, wx.ALIGN_TOP | wx.LEFT, border=20)
        hbox4.Add(self.input_expression_3d, 211, wx.EXPAND | wx.RIGHT, border=20)
        hbox4.Add(label_expression_2d, 110, wx.ALIGN_TOP | wx.LEFT, border=20)
        hbox4.Add(self.input_expression_2d, 210, wx.EXPAND | wx.RIGHT, border=20)

        hbox7 = wx.BoxSizer()
        hbox7.Add(label_curve_point_3d, 113, wx.ALIGN_TOP | wx.LEFT, border=20)
        hbox7.Add(self.input_curve_point_3d, 213, wx.EXPAND | wx.RIGHT, border=20)
        hbox7.Add(label_curve_point_2d, 112, wx.ALIGN_TOP | wx.LEFT, border=20)
        hbox7.Add(self.input_curve_point_2d, 212, wx.EXPAND | wx.RIGHT, border=20)

        hbox5 = wx.BoxSizer()
        hbox5.Add(self.button_APS_3d, 1, flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=5)
        hbox5.Add(self.button_APS_2d, 1, flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=5)

        hbox6 = wx.BoxSizer()
        hbox6.Add(self.button_Clear, 1, flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=5)
        hbox6.Add(self.button_Back, 1, flag=wx.ALIGN_CENTER | wx.LEFT | wx.RIGHT, border=5)

        # 创建垂直方向box布局管理器
        vbox = wx.BoxSizer(wx.VERTICAL)
        # 将hbox添加到vbox
        vbox.Add(hbox1, proportion=1, flag=wx.CENTER | wx.ALL | wx.EXPAND, border=10)
        vbox.Add(hbox2, proportion=1, flag=wx.CENTER | wx.ALL | wx.EXPAND, border=10)
        vbox.Add(hbox3, proportion=1, flag=wx.CENTER | wx.ALL | wx.EXPAND, border=10)
        vbox.Add(hbox4, proportion=1, flag=wx.CENTER | wx.ALL | wx.EXPAND, border=10)
        vbox.Add(hbox7, proportion=1, flag=wx.CENTER | wx.ALL | wx.EXPAND, border=10)
        vbox.Add(hbox5, proportion=1, flag=wx.CENTER | wx.TOP | wx.EXPAND, border=10)
        vbox.Add(hbox6, proportion=1, flag=wx.CENTER | wx.BOTTOM | wx.EXPAND, border=10)
        # 整个界面为一个面板，面板中设置一个垂直方向的布局管理器（根布局管理器）
        panel.SetSizer(vbox)

        # fgs1 = wx.FlexGridSizer(3, 6, 5, 5)
        # fgs1.AddMany([
        #     (label_dna_size, 101, wx.ALIGN_CENTER), (self.input_dna_size, 201, wx.ALIGN_CENTER),
        #     (label_pop_size, 102, wx.ALIGN_CENTER), (self.input_pop_size, 202, wx.ALIGN_CENTER),
        #     (label_accuracy, 109, wx.ALIGN_CENTER), (self.input_accuracy, 209, wx.ALIGN_CENTER),
        #
        #     (label_x_bound, 106, wx.ALIGN_CENTER), (self.input_x_bound, 206, wx.ALIGN_CENTER),
        #     (label_y_bound, 107, wx.ALIGN_CENTER), (self.input_y_bound, 207, wx.ALIGN_CENTER),
        #     (label_z_bound, 108, wx.ALIGN_CENTER), (self.input_z_bound, 208, wx.ALIGN_CENTER),
        #
        #     (label_crossover_rate, 103, wx.ALIGN_CENTER), (self.input_crossover_rate, 203, wx.ALIGN_CENTER),
        #     (label_mutation_rate, 104, wx.ALIGN_CENTER), (self.input_mutation_rate, 204, wx.ALIGN_CENTER),
        #     (label_n_generations, 105, wx.ALIGN_CENTER), (self.input_n_generations, 205, wx.ALIGN_CENTER)
        # ])
        # # 指定行，（行索引，行所占空间比例）
        # fgs1.AddGrowableRow(0, 1)  # 高度占1/13
        # fgs1.AddGrowableRow(1, 1)  # 占1/13
        # fgs1.AddGrowableRow(2, 1)  # 占1/13
        # # 指定列，（列索引，列所占空间比例）
        # fgs1.AddGrowableCol(0, 1)
        # fgs1.AddGrowableCol(1, 2)
        # fgs1.AddGrowableCol(2, 1)
        # fgs1.AddGrowableCol(3, 2)
        # fgs1.AddGrowableCol(4, 1)
        # fgs1.AddGrowableCol(5, 2)
        #
        # fgs2 = wx.FlexGridSizer(1, 4, 5, 5)
        # fgs2.AddMany([
        #     (label_expression_2d, 110, wx.ALIGN_TOP), (self.input_expression_2d, 210, wx.EXPAND),
        #     (label_expression_3d, 111, wx.ALIGN_TOP), (self.input_expression_3d, 211, wx.EXPAND)
        # ])
        # # 指定行，（行索引，行所占空间比例）
        # fgs2.AddGrowableRow(0, 1)
        # # 指定列，（列索引，列所占空间比例）
        # fgs2.AddGrowableCol(0, 1)
        # fgs2.AddGrowableCol(1, 2)
        # fgs2.AddGrowableCol(2, 1)
        # fgs2.AddGrowableCol(3, 2)
        #
        # fgs3 = wx.FlexGridSizer(1, 4, 5, 5)
        # fgs3.AddMany([
        #     (self.button_APS_3d, 1, wx.ALIGN_CENTER), (self.button_APS_2d, 2, wx.ALIGN_CENTER),
        #     (self.button_Clear, 3, wx.ALIGN_CENTER), (self.button_Back, 4, wx.ALIGN_CENTER)
        # ])
        # # 指定行，（行索引，行所占空间比例）
        # fgs3.AddGrowableRow(0, 1)
        # # 指定列，（列索引，列所占空间比例）
        # fgs3.AddGrowableCol(0, 1)
        # fgs3.AddGrowableCol(1, 1)
        # fgs3.AddGrowableCol(2, 1)
        # fgs3.AddGrowableCol(3, 1)
        #
        # # 将FlexGrid布局管理器添加进面板
        # vbox = wx.BoxSizer(wx.VERTICAL)
        # vbox.Add(fgs1, proportion=1, flag=wx.CENTER | wx.EXPAND, border=10)
        # vbox.Add(fgs2, proportion=1, flag=wx.CENTER | wx.EXPAND, border=10)
        # # vbox.Add(self.label_result, proportion=1, flag=wx.CENTER, border=5)
        # vbox.Add(fgs3, proportion=1, flag=wx.CENTER | wx.EXPAND, border=10)
        # panel.SetSizer(vbox)

    def OnButton_APS_3d(self, event):
        APS_3d.DNA_SIZE = int(self.input_dna_size.GetValue())
        APS_3d.POP_SIZE = int(self.input_pop_size.GetValue())
        APS_3d.CROSSOVER_RATE = float(self.input_crossover_rate.GetValue())
        APS_3d.MUTATION_RATE = float(self.input_mutation_rate.GetValue())
        APS_3d.N_GENERATIONS = int(self.input_n_generations.GetValue())

        pattern = re.compile(r'-*\d+')
        listX = pattern.findall(self.input_x_bound.GetValue())
        APS_3d.X_BOUND = [int(i) for i in listX]
        APS_3d.X_BOUND.sort()
        listY = pattern.findall(self.input_y_bound.GetValue())
        APS_3d.Y_BOUND = [int(j) for j in listY]
        APS_3d.Y_BOUND.sort()
        listZ = pattern.findall(self.input_z_bound.GetValue())
        APS_3d.Z_BOUND = [int(k) for k in listZ]
        APS_3d.Z_BOUND.sort()

        APS_3d.ACCURACY = int(self.input_accuracy.GetValue())
        APS_3d.EXPRESSION_3D = self.input_expression_3d.GetValue()
        APS_3d.CURVE_POINT_3D = self.input_curve_point_3d.GetValue()

        APS_3d.print_3d()
        event.Skip()

    def OnButton_APS_2d(self, event):
        APS_2d.DNA_SIZE = int(self.input_dna_size.GetValue())
        APS_2d.POP_SIZE = int(self.input_pop_size.GetValue())
        APS_2d.CROSSOVER_RATE = float(self.input_crossover_rate.GetValue())
        APS_2d.MUTATION_RATE = float(self.input_mutation_rate.GetValue())
        APS_2d.N_GENERATIONS = int(self.input_n_generations.GetValue())

        pattern = re.compile(r'-*\d+')
        listX = pattern.findall(self.input_x_bound.GetValue())
        APS_2d.X_BOUND = [int(i) for i in listX]
        APS_2d.X_BOUND.sort()
        listY = pattern.findall(self.input_y_bound.GetValue())
        APS_2d.Y_BOUND = [int(j) for j in listY]
        APS_2d.Y_BOUND.sort()

        APS_2d.ACCURACY = int(self.input_accuracy.GetValue())
        APS_2d.EXPRESSION_2D = self.input_expression_2d.GetValue()
        APS_2d.CURVE_POINT_2D = self.input_curve_point_2d.GetValue()

        APS_2d.print_2d()
        event.Skip()

    def OnButton_Clear(self, event):
        self.input_dna_size.SetValue(str(DNA_SIZE))
        self.input_pop_size.SetValue(str(POP_SIZE))
        self.input_crossover_rate.SetValue(str(CROSSOVER_RATE))
        self.input_mutation_rate.SetValue(str(MUTATION_RATE))
        self.input_n_generations.SetValue(str(N_GENERATIONS))

        self.input_x_bound.SetValue(str(X_BOUND))
        self.input_y_bound.SetValue(str(Y_BOUND))
        self.input_z_bound.SetValue(str(Z_BOUND))

        self.input_accuracy.SetValue(str(ACCURACY))
        self.input_expression_2d.SetValue(str(EXPRESSION_2D))
        self.input_expression_3d.SetValue(str(EXPRESSION_3D))
        self.input_curve_point_2d.SetValue(str(CURVE_POINT_2D))
        self.input_curve_point_3d.SetValue(str(CURVE_POINT_3D))
        event.Skip()

    def OnButton_Back(self, event):
        dlg = wx.MessageDialog(self, '确认返回？', '返回', wx.CANCEL | wx.OK | wx.ICON_QUESTION)
        if dlg.ShowModal() == wx.ID_OK:
            self.Destroy()


# 自定义应用程序对象
class App(wx.App):
    def OnInit(self):
        # 创建窗口对象
        frameAPS = APSFrame(None)
        frameAPS.Show()
        return True

    # def OnExit(self):
    #     print('应用程序退出')
    #     self.Destroy()
    #     return 0


if __name__ == '__main__':
    app = App()  # 调用上面函数
    app.MainLoop()  # 进入主事件循环
