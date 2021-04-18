import wx
import GA_GUI
import APS_GUI
import sys


# 自定义窗口类MyFrame
class MainFrame(wx.Frame):
    def __init__(self, superior):
        wx.Frame.__init__(self, parent=superior, title=u'GUI界面')
        self.Centre()
        self.SetSize(300, 150)
        panel = wx.Panel(self, -1)

        # 添加按钮
        self.buttonGA = wx.Button(parent=panel, label='遗传算法')
        self.buttonAPS = wx.Button(parent=panel, label='APS排程')
        self.buttonQuit = wx.Button(parent=panel, label='退出')

        # 绑定事件处理函数
        self.Bind(wx.EVT_BUTTON, handler=self.OnButtonGA, source=self.buttonGA)
        self.Bind(wx.EVT_BUTTON, handler=self.OnButtonAPS, source=self.buttonAPS)
        self.Bind(wx.EVT_BUTTON, handler=self.OnButtonQuit, source=self.buttonQuit)

        self.panel = panel
        fgs1 = wx.FlexGridSizer(1, 2, 10, 10)
        fgs1.Add(self.buttonGA, wx.ALIGN_CENTER)
        fgs1.Add(self.buttonAPS, wx.ALIGN_CENTER)

        fgs2 = wx.FlexGridSizer(1, 1, 10, 10)
        fgs2.Add(self.buttonQuit)

        # 将FlexGrid布局管理器添加进面板
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(fgs1, proportion=2, flag=wx.CENTER | wx.TOP, border=20)
        vbox.Add(fgs2, proportion=1, flag=wx.CENTER | wx.BOTTOM, border=20)
        panel.SetSizer(vbox)

    def OnButtonGA(self, event):
        app_GA = GA_GUI.App()
        app_GA.MainLoop()

    def OnButtonAPS(self, event):
        app_APS = APS_GUI.App()
        app_APS.MainLoop()

    def OnButtonQuit(self, event):
        dlg = wx.MessageDialog(self, '确认退出？', '退出', wx.CANCEL | wx.OK | wx.ICON_QUESTION)
        if dlg.ShowModal() == wx.ID_OK:
            print('应用程序退出')
            self.Destroy()
            sys.exit(0)


# 自定义应用程序对象
class App(wx.App):
    def OnInit(self):
        # 创建窗口对象
        framemain = MainFrame(None)
        framemain.Show()
        return True

    def OnExit(self):
        print('应用程序退出')
        sys.exit(0)
        return 0


if __name__ == '__main__':
    app = App()  # 调用上面函数
    app.MainLoop()  # 进入主事件循环
