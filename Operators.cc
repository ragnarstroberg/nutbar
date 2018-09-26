
#include "Operators.hh"
#include "Settings.hh"







ScalarNME::ScalarNME(Settings& settings)
{
  OneBody.set_size(settings.J2_f.size(), settings.J2_i.size() );
  TwoBody.set_size(settings.J2_f.size(), settings.J2_i.size() );
  for (size_t indexJi=0;indexJi<settings.J2_i.size();++indexJi )
  {
   for (size_t indexJf=0;indexJf<settings.J2_f.size();++indexJf )
   {
     if (settings.J2_f[indexJf]==settings.J2_i[indexJi])
     {
        OneBody(indexJf,indexJi).zeros(settings.NJ_f[indexJf],settings.NJ_i[indexJi]);
        TwoBody(indexJf,indexJi).zeros(settings.NJ_f[indexJf],settings.NJ_i[indexJi]);
     }
   }
  }
}






TensorNME::TensorNME(Settings& settings, TensorOperator& TensorOp)
 : Lambda(TensorOp.Lambda)
{
  OneBody.set_size(settings.J2_f.size(), settings.J2_i.size() );
  TwoBody.set_size(settings.J2_f.size(), settings.J2_i.size() );
  for (size_t indexJi=0;indexJi<settings.J2_i.size();++indexJi )
  {
   for (size_t indexJf=0;indexJf<settings.J2_f.size();++indexJf )
   {
//     if ( std::abs(settings.J2_i[indexJi] - settings.J2_f[indexJf])> 2*TensorOp.Lambda) continue;

     OneBody(indexJf,indexJi).zeros(settings.NJ_f[indexJf],settings.NJ_i[indexJi]);
     TwoBody(indexJf,indexJi).zeros(settings.NJ_f[indexJf],settings.NJ_i[indexJi]);
   }
  }

}




DaggerNME::DaggerNME(Settings& settings, DaggerOperator& DaggerOp)
 : Lambda2(DaggerOp.Lambda2)
{
  ax.set_size(settings.J2_f.size(), settings.J2_i.size() );
  axaxa.set_size(settings.J2_f.size(), settings.J2_i.size() );
  for (size_t indexJi=0;indexJi<settings.J2_i.size();++indexJi )
  {
   for (size_t indexJf=0;indexJf<settings.J2_f.size();++indexJf )
   {
//     if ( std::abs(settings.J2_i[indexJi] - settings.J2_f[indexJf])> DaggerOp.Lambda2) continue;

     ax(indexJf,indexJi).zeros(settings.NJ_f[indexJf],settings.NJ_i[indexJi]);
     axaxa(indexJf,indexJi).zeros(settings.NJ_f[indexJf],settings.NJ_i[indexJi]);
   }
  }

}
