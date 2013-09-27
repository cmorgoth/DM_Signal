 {
   gROOT->Reset();
   std::ifstream mfile0("list_of_files_table.list");
   std::ofstream outfile("table_eff.tex");
   outfile << "\\begin{table}[htdp]\n\\caption{default}\n\\begin{center}\n\\begin{tabular}{|c|c|c|c|c|c|c|}\n\\hline\n";
   outfile << "Sample" << " & " << "$R^{2} > 0.5$" << " & " << "$R^{2} > 0.5 && \\Delta\\Phi < 2.5$" << " & " << "$MR*R^{2} > 100$" << " & " << "$MR*R^{2} > 100$" << " & " << "$MR*R^{2} > 100$" << " & " << "$MR*R^{2} > 100$" << "\\" << "\\" << "\n";
     outfile << "\\hline" << std::endl;
   std::string fname0;
   std::cout.precision(16);
   int ctr = 0, ctr2 =0;
   double arr[7][27];
   std::string ss[27];
   if (mfile0.is_open()){
     while ( mfile0.good() ){
       mfile0 >> fname0;
       std::cout << fname0 << std::endl;
       std::ifstream mfile1(fname0.c_str());
       ctr = 0; 
       if (mfile1.is_open()){
	 while ( mfile1.good() ){
	   mfile1 >> ss[ctr] >> arr[ctr2][ctr];
	   if(mfile1.eof())break;
	   ctr++;
	   std::cout << "name: " << ss[ctr] << " " << arr[ctr2][ctr] << " " << ctr << " " << ctr2<< std::endl;
	 }
       }else{
	 std::cout << "unable to open the file" << std::endl;
       }
       mfile1.close();
       if(mfile0.eof())break;
       ctr2++;
     }
   }else{
     std::cout << "Unable to opne the file" << std::end; 
   }
   
 
   for(int j = 0; j < 24; j++){
     outfile << ss[j] << " & " << arr[0][j] << " & " << arr[1][j] << " & " << arr[2][j] << " & " << arr[3][j] << " & " << arr[4][j] << " & " << arr[5][j] << "\\%" << "\\" << "\\" << "\n";
     outfile << "\\hline" << std::endl;
     
    }
   outfile << "\\end{tabular}\n\\end{center}\n\\label{default}\n\\end{table}\n";
   outfile.close();
}     
