#include "interprete.hh"

interprete::interprete()
{
  eval=true;
  cutoff=1.0e-6;
}

void
interprete::set_cutoff(double r)
{
  cutoff=r;
}

double
interprete::S_to_Val(string *str)
{
  double v=atof((*str).c_str());
#ifdef PRINT_IT
  cout << "v=" << v << "\n";
#endif
  return(v);
}

double
interprete::get_value(string *str, double ll)
{
#ifdef PRINT_IT
  cout << "value[" << *str << "]=" << variable[(*str)].value << "\n";
#endif
  return(variable[(*str)].value);
}

void
interprete::put_value(string *str, int num,int Type, double val)
{
#ifdef PRINT_IT
  cout << "*str=" << *str << " val=" << val << "\n";
#endif
  tmp_put_value_name=*str;
  variable[(*str)].value=val;
  variable[(*str)].index=num;
  variable[(*str)].type=Type;
}

void
interprete::print_all()
{
  map<string, t_val_index>::iterator iter;
  for( iter = variable.begin(); iter != variable.end(); iter++ )
    {
      cout << "Name= " << iter->first << ", Value= " << iter->second.value << ", Index= " <<  iter->second.index << ", Type= " <<  iter->second.type << ", Indexed=" << iter->second.index << endl;
    }
}

void
interprete::create_id_map(int *Table, int Size, int HSize)
{
  map<string, t_val_index>::iterator iter;
  for( iter = variable.begin(); iter != variable.end(); iter++ )
    {
      if(iter->second.type==0)
        {
          i_endo_variable[iter->second.index].name=iter->first;
          i_endo_variable[iter->second.index].value=iter->second.value;
        }
      else if(iter->second.type==1)
        {
          i_exo_variable[iter->second.index].name=iter->first;
          i_exo_variable[iter->second.index].value=iter->second.value;
        }
      else if(iter->second.type==4)
        {
          i_param_variable[iter->second.index].name=iter->first;
          i_param_variable[iter->second.index].value=iter->second.value;
        }
    }
}

double
interprete::GetDataValue(int id, Type type)
{
  switch(type)
    {
    case eParameter:
      return(i_param_variable[(long int) id].value);
    case eEndogenous:
      return(i_endo_variable[(long int) id].value);
    case eExogenous:
      return(i_exo_variable[(long int) id].value);
    default:
      cerr << "interprete::GetDataValue: unhandled type!!";
      exit(-1);
    }
}
