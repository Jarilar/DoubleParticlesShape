using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using MahApps.Metro.Controls;

namespace DoubleP
{
    /// <summary>
    /// MainWindow.xaml 的交互逻辑
    /// </summary>
    public partial class MainWindow : MetroWindow
    {
        public MainWindow()
        {
            InitializeComponent();
        }
        public int num;//行数计数
        public int Nocol;//当前所在行数
        public int DRHdata,RoriginM;
        public double RLdry, RBdry,Ria,Rib,kdata,BRi,Vs;
        public double[,] RMB1, RML1,RMB2;
        public double[] RL, RB, Rgd, Vgd,RML2,RML3;
        
        public double[] NSO = new double[]{1.4510,1.4620,1.4830,1.5060,1.5310,1.5580,1.5870,
                                             1.6200,1.6560,1.6960,1.7410,1.7930,1.8520,1.9200,
                                              2.0020,2.1010,2.2260,2.39,2.6220,2.9910,3.7510};
        public double[] NCL = new double[]{1.6340,1.6570,1.6810 ,1.7070 ,1.7340 ,1.7630 ,1.7940,
                                         1.8270, 1.8630, 1.9020 ,1.9430 ,1.9890 ,2.0400 ,2.0960,
                                        2.1580, 2.2280 ,2.3090 ,2.4020 ,2.5120 ,2.6450 ,2.8110,
                                        3.0290 ,3.3330 ,3.8150 ,4.8020};
        public double[] NaSO = new double[] {1.6280,1.6590,1.6930,1.7290,1.7700,1.8150,1.8660,
                                              1.9230,1.9880,2.0650,2.1550,2.2650,2.4030,2.5840,
                                              2.8380,3.2420,4.0710};
        public double[] NNO = new double[] {1.4360,1.4510,1.4680,1.4850,1.5040,1.5230,1.5440,
                                             1.5650,1.5890,1.6130,1.6400,1.6680,1.6990,1.7320,
                                             1.7680,1.8080,1.8520,1.9000,1.9540,2.0150,2.0860,
                                             2.1670,2.2640,2.3810,2.5280,2.7200,2.9900,3.4180,4.2960 };


        private void MetroWindow_Closing(object sender, System.ComponentModel.CancelEventArgs e)
        {
            MessageBoxResult result = MessageBox.Show("确认是否要退出？", "退出", MessageBoxButton.YesNo);//显示确认窗口
            if (result == MessageBoxResult.No)
            {
                e.Cancel = true;//取消操作
            }
        }

        private void MetroWindow_Closed(object sender, EventArgs e)
        {
            Application.Current.Shutdown();//先停止线程,然后终止进程.
            Environment.Exit(0);//直接终止进程.
        }


        #region Buttons

        /// <summary>
        /// 潮解前
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void BeShapeBtn_Click(object sender, RoutedEventArgs e)
        {

        }

        private void BeDielBtn_Click(object sender, RoutedEventArgs e)
        {

        }

        private void BeParfile_Click(object sender, RoutedEventArgs e)
        {

        }
        /// <summary>
        /// 潮解后未完全没入
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void OutShapeBtn_Click(object sender, RoutedEventArgs e)
        {

        }

        private void OutDielBtn_Click(object sender, RoutedEventArgs e)
        {

        }

        private void OutParfile_Click(object sender, RoutedEventArgs e)
        {

        }

        /// <summary>
        /// 潮解后完全没入
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void InShapeBtn_Click(object sender, RoutedEventArgs e)
        {

        }

        private void InDielBtn_Click(object sender, RoutedEventArgs e)
        {

        }

        private void InParfile_Click(object sender, RoutedEventArgs e)
        {

        }
        /// <summary>
        /// 亲水性数据录入
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void D1Btn_Click(object sender, RoutedEventArgs e)
        {

            kdata = double.Parse(K_Data.Text);
            DRHdata = int.Parse(DRH.Text);
            RLdry = double.Parse(RLorgin.Text);
            RoriginM = int.Parse(R1orgin.Text);

            string[] Datams = Ri1.Text.Split('E');
            Ria = double.Parse(Datams[0]);
            Rib = double.Parse(Datams[1]);


            LAHM(kdata, DRHdata, RLdry);



            Datashow.Document.Blocks.Add(new Paragraph(new Run("亲水性气溶胶粒子运算完成！")));
            switch (kdata)
            {
                case 1.12:
                    Datashow.Document.Blocks.Add(new Paragraph(new Run("物质为：氯化钠(Sodium chloride)")));
                    break;
                case 0.68:
                    Datashow.Document.Blocks.Add(new Paragraph(new Run("物质为：硫酸钠(Sodium sulfate)")));
                    break;
                case 0.8:
                    Datashow.Document.Blocks.Add(new Paragraph(new Run("物质为：硝酸钠(Sodium nitrate)")));
                    break;
                case 0.53:
                    Datashow.Document.Blocks.Add(new Paragraph(new Run("物质为：硫酸铵(Ammonium sulfate)")));
                    break;
            }
            Datashow.Document.Blocks.Add(new Paragraph(new Run("干粒径为："+RLorgin.Text)));
            Datashow.Document.Blocks.Add(new Paragraph(new Run("折射率变化方程：y=" +Ria.ToString()+"x+"+Rib.ToString())));

        }

        /// <summary>
        /// 疏水性数据录入
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void D2Btn_Click(object sender, RoutedEventArgs e)
        {
            RBdry = double.Parse(RBorgin.Text) ;           
            BRi= double.Parse(Ri2.Text);

            BAHM(RBdry,DRHdata);

            Datashow.Document.Blocks.Add(new Paragraph(new Run("疏水性气溶胶粒子运算完成！")));
            Datashow.Document.Blocks.Add(new Paragraph(new Run("预祝计算顺利！")));
        }

        /// <summary>
        /// 亲水性数据清除
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void Cl1Btn_Click(object sender, RoutedEventArgs e)
        {
            K_Data.Text = "";
            DRH.Text = "";
            R1orgin.Text = "";
            RLorgin.Text = "";
            Ri1.Text = "";
        }
        /// <summary> 
        /// 疏水性数据清除
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void Cl2Btn_Click(object sender, RoutedEventArgs e)
        {
            R2orgin.Text = "";
            RBorgin.Text = "";
            Ri2.Text = "";
        }


        #endregion

        #region Functions
        /// <summary>
        /// LAHM算法
        /// </summary>
        /// <param name="k" 吸湿性参数></param>
        /// <param name="DRH" 潮解相对湿度></param>
        /// <param name="R" 干粒径></param>
        private void LAHM(double k, int DRH, double R)
        {
            int Arrsize = 100 - DRH;
            double[] GRH = new double[Arrsize];
            RL = new double[101];
            double Vdry = Math.PI * 4 * Math.Pow(R, 2) / 3;
            for (int i = 0; i <= DRH; i++)
            {
                RL[i] = 3 * Math.Pow(((Vdry * k * i / (100 - i)) + Vdry), 1 / 3) / (4 * Math.PI);
            }
            switch (k)
            {
                case 1.12:
                    GRH = NCL;
                    break;
                case 0.68:
                    GRH = NaSO;
                    break;
                case 0.8:
                    GRH = NNO;
                    break;
                case 0.53:
                    GRH = NSO;
                    break;
            }
            for (int j = 0; j < Arrsize; j++)
            {
                RL[j + DRH + 1] = GRH[j] * R;
            }
        }
        /// <summary>
        /// BAHM算法
        /// </summary>
        /// <param name="R" 干粒径></param>
        /// <param name="DRH" 潮解相对湿度></param>
        private void BAHM(double R, int DRH)
        {
            RB = new double[101];
            for (int i = 1; i <= DRH; i++)
            {
                double m = Math.Pow(4.34 / (Math.Log(100 / i) + (1.08 / R)), 1 / 3) * 0.00004;
                RB[i] = m + R;
            }
            for (int j = DRH; j <100; j++)
            {
                double m = Math.Pow(4.34 / (Math.Log(100 / j) + (1.08 / R)), 1 / 3) * 0.00004;
                RB[j+1] = m + R;
            }
        }
        
        /// <summary>
        /// 潮解前模型半径
        /// </summary>
        /// <param name="DRH"></param>
        private void BeModel (int DRH)
        {
            int Rm0 = RoriginM;
            Vs = Rm0 / RLdry;
            RML1 = new double[DRH, 2];
            RMB1 = new double[DRH, 2];

            for(int i=0;i<=DRH;i++)
            {
                RML1[i, 0] = Rm0;
                RML1[i, 1] = RL[i] * Vs;
                RMB1[i, 0] = Rm0;
                RMB1[i, 1] = RB[i] * Vs;
            }
        }

        /// <summary>
        /// 半吞没形态模型半径
        /// </summary>
        /// <param name="DRH"></param>
        /// <param name="Distan" 实际球心距></param>
        private void OutModel(int DRH,double Distan)
        {
            int Arrsize = 100 - DRH;
            RML2 = new double[Arrsize];
            RMB2 = new double[Arrsize, 2];
            Vgd = new double[Arrsize];
            double Jdian = ((Distan * Distan) - (RBdry * RBdry) + (RL[0] * RL[0])) / (2 * Distan);
            double h2 = Jdian - (Distan - RBdry);

            for (int i = 0; i < Arrsize; i++)
            {
                RMB2[i, 0] = RoriginM;
                RMB2[i, 1] = RB[i] * Vs;
                double h1 = RL[i] - Jdian;
                Vgd[i] = ((Math.PI * h1 * (3 * Rgd[i] * Rgd[i] + h1 * h1)) / 6) + ((Math.PI * h2 * (3 * RBdry * RBdry + h2 * h2)) / 6) + (RL[i] * RL[i] * RL[i] * 4 * Math.PI / 3);
                Rgd[i]=
                RML2[i] = (Math.Pow(Vgd[i], 1 / 3) * 3 / (4 * Math.PI)) *Vs;
            }
        }

        private void InsieModel(int DRH,int RMorigin)
        {

        }

        private void Demition(double Distence, int Nnumber, double R00)
        {
            /*
            double Rmax, DRmax;
            int n = 1;//计数用
            string[,] arr = new string[500000, 7];//定义一个字符数组

            ////录入数据////////////////////////////////////////////

            //获得最大半径           

            //生成形状

            arr[0, 0] = ">TARCEL:";
            arr[0, 1] = "shape;";
            n++;

            arr[1, 1] = "=";
            arr[1, 2] = "NAT";
            n++;
            arr[2, 1] = "1.000000";
            arr[2, 2] = "0.000000";
            arr[2, 3] = "0.000000";
            arr[2, 4] = "=";
            arr[2, 5] = "A_1 vector";
            n++;
            arr[3, 0] = "0.000000";
            arr[3, 1] = "1.000000";
            arr[3, 2] = "0.000000";
            arr[3, 3] = "=";
            arr[3, 4] = "A_2 vector";
            n++;
            arr[4, 0] = "1.000000";
            arr[4, 1] = "1.000000";
            arr[4, 2] = "1.000000";
            arr[4, 3] = "=";
            arr[4, 4] = "lattice spacings (d_x,d_y,d_z)/d";
            n++;
            arr[5, 0] = "-0.500000";
            arr[5, 1] = "-0.500000";
            arr[5, 2] = "-0.500000";
            arr[5, 3] = "=";
            arr[5, 4] = "lattice offset x0(1-3) ";
            arr[5, 5] = "=";
            arr[5, 6] = "(x_TF,y_TF,z_TF)/d for dipole 0 0 0";
            n++;
            arr[6, 0] = "JA";
            arr[6, 1] = "IX";
            arr[6, 2] = "IY";
            arr[6, 3] = "IZ";
            arr[6, 4] = "ICOMP(x,y,z)";

            double R1 = RML[Nnumber];
            double R2 = RMB[Nnumber, 0];
            double R3 = RMB[Nnumber, 1];

            Rmax = R3 + R1;

            DRmax = 2 * Rmax;
            for (int i = -(int)DRmax; i <= DRmax; i++)
            {
                for (int j = -(int)Rmax; j <= Rmax; j++)
                {
                    for (int k = -(int)Rmax; k <= Rmax; k++)
                    {
                        double CA = Math.Sqrt(Math.Pow(i, 2) + Math.Pow(j, 2) + Math.Pow(k, 2));
                        double CB = Math.Sqrt(Math.Pow(i - R00, 2) + Math.Pow(j, 2) + Math.Pow(k, 2));

                        if (CA <= R1 && CB > R2)//亲水性的
                        {
                            arr[n, 0] = (n - 6).ToString();
                            arr[n, 1] = i.ToString();
                            arr[n, 2] = j.ToString();
                            arr[n, 3] = k.ToString();
                            arr[n, 4] = "3";
                            arr[n, 5] = "3";
                            arr[n, 6] = "3";
                            n++;
                        }
                        else if (CB <= R2)//疏水性的
                        {
                            arr[n, 0] = (n - 7).ToString();
                            arr[n, 1] = i.ToString();
                            arr[n, 2] = j.ToString();
                            arr[n, 3] = k.ToString();
                            arr[n, 4] = "2";
                            arr[n, 5] = "2";
                            arr[n, 6] = "2";
                            n++;
                        }
                        else if (CB > R2 && CB <= R3 && CA > R1)//表面水层的
                        {
                            arr[n, 0] = (n - 7).ToString();
                            arr[n, 1] = i.ToString();
                            arr[n, 2] = j.ToString();
                            arr[n, 3] = k.ToString();
                            arr[n, 4] = "1";
                            arr[n, 5] = "1";
                            arr[n, 6] = "1";
                            n++;
                        }
                    }
                }

                arr[1, 0] = (n - 7).ToString();
                //写入txt文档
                string Fname = Nnumber.ToString() + "_" + Distence.ToString();
                string ModelName = "双球体";
                ShowBox.Document.Blocks.Add(new Paragraph(new Run("模型：" + ModelName + "RH=" + Fname + "%     完成" + "\r\n")));
                //保存地址为D盘自定义文件名
                FileStream fs = new FileStream("D:\\test\\" + Fname + ".txt", FileMode.Create);
                StreamWriter sw = new StreamWriter(fs);
                //初始化二维数组s 用来接收arr
                string[,] s = new string[500000, 7];
                //接收arr
                s = arr;

                for (int l = 0; l < 500000; ++l)
                {
                    for (int h = 0; h < 7; ++h)
                    {
                        //s的每个值和arr的每个值对应
                        s[l, h] = arr[l, h];

                        string output;

                        output = s[l, h];
                        //有个空格作为间隔，
                        sw.Write(output + " ");
                    }
                    sw.WriteLine();
                }
                //清空缓冲区
                sw.Flush();
                //关闭流
                sw.Close();
                fs.Close();
            }
            */
        }
        //生成配置文件
        private void writeparfile(double distand, int numb, double aeff)
        {
            /*
            int n = 1;//计数用
            string[,] arr = new string[500000, 7];//定义一个字符数组
            arr[0, 0] = "' ========== Parameter file for v7.3 ==================='";
            arr[1, 0] = "'**** Preliminaries ****' ";
            arr[2, 0] = "'NOTORQ' = CMDTRQ*6 (DOTORQ, NOTORQ) -- either do or skip torque calculations";
            arr[3, 0] = "'PBCGS2' = CMDSOL*6 (PBCGS2, PBCGST, GPBICG, QMRCCG, PETRKP) -- CCG method ";
            arr[4, 0] = "'GPFAFT' = CMETHD*6 (GPFAFT, FFTMKL) -- FFT method";
            arr[5, 0] = "'GKDLDR' = CALPHA*6 (GKDLDR, LATTDR, FLTRCD) -- DDA method ";
            arr[6, 0] = "'NOTBIN' = CBINFLAG (NOTBIN, ORIBIN, ALLBIN) -- binary output? ";
            arr[7, 0] = "'**** Initial Memory Allocation ****' ";
            arr[8, 0] = "200 200 200 = dimensioning allowance for target generation ";
            arr[9, 0] = "'**** Target Geometry and Composition ****' ";
            arr[10, 0] = "'FRMFILPBC' = CSHAPE*9 shape directive";
            arr[11, 0] = "0 0 '../shapefile/" + numb.ToString() + "_" + distand.ToString() + ".txt'";
            arr[12, 0] = "2         = NCOMP = number of dielectric materials";
            arr[13, 0] = "'../diel/" + PlanArray[1, 4].ToString() + ".txt' = file with refractive index 1";
            arr[14, 0] = "'../diel/" + PlanArray[1, 5].ToString() + ".txt' = file with refractive index 2";
            arr[15, 0] = "'**** Additional Nearfield calculation? ****' ";
            arr[16, 0] = "0 = NRFLD (=0 to skip nearfield calc., =1 to calculate nearfield E) ";
            arr[17, 0] = "0.0 0.0 0.0 0.0 0.0 0.0 (fract. extens. of calc. vol. in -x,+x,-y,+y,-z,+z)";
            arr[18, 0] = "'**** Error Tolerance ****'";
            arr[19, 0] = "1.00e-5 = TOL = MAX ALLOWED (NORM OF |G>=AC|E>-ACA|X>)/(NORM OF AC|E>)";
            arr[20, 0] = "'**** Maximum number of iterations ****'";
            arr[21, 0] = "10000     = MXITER";
            arr[22, 0] = "'**** Integration limiter for PBC calculations ****'";
            arr[23, 0] = "1.00e-2 = GAMMA (1e-2 is normal, 3e-3 for greater accuracy)";
            arr[24, 0] = "'**** Angular resolution for calculation of <cos>, etc. ****'";
            arr[25, 0] = "0.5	= ETASCA (number of angles is proportional to [(3+x)/ETASCA]^2 )";
            arr[26, 0] = "'**** Wavelengths (micron) ****'";
            arr[27, 0] = "0.55 0.55 1 'INV' = wavelengths (1st,last,howmany,how=LIN,INV,LOG,TAB)";
            arr[28, 0] = "'**** Refractive index of ambient medium ****'";
            arr[29, 0] = "1.0000 = NAMBIENT";
            arr[30, 0] = "'**** Effective Radii (micron) **** '";
            arr[31, 0] = aeff.ToString();
            arr[31, 1] = aeff.ToString() + " 1 'LIN' = eff. radii (1st,last,howmany,how=LIN,INV,LOG,TAB)";
            arr[32, 0] = "'**** Define Incident Polarizations ****'";
            arr[33, 0] = "(0,0) (1.,0.) (0.,0.) = Polarization state e01 (k along x axis)";
            arr[34, 0] = "2 = IORTH  (=1 to do only pol. state e01; =2 to also do orth. pol. state)";
            arr[35, 0] = "'**** Specify which output files to write ****'";
            arr[36, 0] = "1 = IWRKSC (=0 to suppress, =1 to write .sca file for each target orient.";
            arr[37, 0] = "'**** Specify Target Rotations ****'";
            arr[38, 0] = "0.    0.   1  = BETAMI, BETAMX, NBETA  (beta=rotation around a1)";
            arr[39, 0] = "0.    0.   1  = THETMI, THETMX, NTHETA (theta=angle between a1 and k)";
            arr[40, 0] = "0.    0.   1  = PHIMIN, PHIMAX, NPHI (phi=rotation angle of a1 around k)";
            arr[41, 0] = "'**** Specify first IWAV, IRAD, IORI (normally 0 0 0) ****'";
            arr[42, 0] = "0   0   0    = first IWAV, first IRAD, first IORI (0 0 0 to begin fresh)";
            arr[43, 0] = "'**** Select Elements of S_ij Matrix to Print ****'";
            arr[44, 0] = "9	= NSMELTS = number of elements of S_ij to print (not more than 9)";
            arr[45, 0] = "11 12 21 22 31 33 44 34 43	= indices ij of elements to print";
            arr[46, 0] = "'**** Specify Scattered Directions ****'";
            arr[47, 0] = "'LFRAME' = CMDFRM (LFRAME, TFRAME for Lab Frame or Target Frame)";
            arr[48, 0] = "1 = NPLANES = number of scattering planes";
            arr[49, 0] = "0.  0. 180. 1 = phi, theta_min, theta_max (deg) for plane A";


            //写入txt文档
            string Fname = "ddcast.par";
            string directoryPath = @"D:\test\" + numb.ToString() + "_" + distand.ToString();//定义一个路径变量
            if (!Directory.Exists(directoryPath))//如果路径不存在
            {
                Directory.CreateDirectory(directoryPath);//创建一个路径的文件夹
            }
            ShowBox.Document.Blocks.Add(new Paragraph(new Run("参数文件：RH=" + Fname + "%     完成" + "\r\n")));
            StreamWriter sw = new StreamWriter(System.IO.Path.Combine(directoryPath, Fname));
            //初始化二维数组s 用来接收arr
            string[,] s = new string[52, 2];
            //接收arr
            s = arr;

            for (int l = 0; l < 52; ++l)
            {
                for (int h = 0; h < 2; ++h)
                {
                    //s的每个值和arr的每个值对应
                    s[l, h] = arr[l, h];

                    string output;

                    output = s[l, h];
                    //有个空格作为间隔，
                    sw.Write(output + " ");
                }
                sw.WriteLine();
            }
            //清空缓冲区
            sw.Flush();
            //关闭流
            sw.Close();
            //fs.Close();
            */
        }
        //生成折射率文件
        private void writediel()
        {
            /*
            string[,] arr = new string[500000, 7];//定义一个字符数组
            for (int i = 4; i < 6; i++)
            {

                string Fname = PlanArray[1, i].ToString() + ".txt";
                arr[0, 0] = "m=" + PlanArray[1, i].ToString();
                arr[1, 0] = "1 2 3 0 0 = columns for wave, Re(n), Im(n), eps1, eps2";
                arr[2, 0] = "   LAMBDA     Re(N)     Im(N)";
                arr[3, 0] = "  0.000001  " + PlanArray[1, i].ToString() + "  0.00000";
                arr[4, 0] = "  1.000000  " + PlanArray[1, i].ToString() + "  0.00000";
                arr[5, 0] = "  100000.0  " + PlanArray[1, i].ToString() + "  0.00000";

                string directoryPath = @"D:\test\diel\";//定义一个路径变量
                if (!Directory.Exists(directoryPath))//如果路径不存在
                {
                    Directory.CreateDirectory(directoryPath);//创建一个路径的文件夹
                }
                monitorbox.Text += "参数文件：RH=" + Fname + "%     完成" + "\r\n";
                StreamWriter sw = new StreamWriter(Path.Combine(directoryPath, Fname));
                //初始化二维数组s 用来接收arr
                string[,] s = new string[52, 2];
                //接收arr
                s = arr;

                for (int l = 0; l < 52; ++l)
                {
                    for (int h = 0; h < 2; ++h)
                    {
                        //s的每个值和arr的每个值对应
                        s[l, h] = arr[l, h];

                        string output;

                        output = s[l, h];
                        //有个空格作为间隔，
                        sw.Write(output + " ");
                    }
                    sw.WriteLine();
                }
                //清空缓冲区
                sw.Flush();
                //关闭流
                sw.Close();
            }
                  */
        }
          

        #endregion
    }
}