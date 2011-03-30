#include "exception.h"
#include <iostream>
#include <sstream> 

namespace exc
{
  
  exception::exception(int code,
		       const char* file,
		       const char* line,
		       const char* reason)
    : m_code(code), m_file(file), m_line(line), m_reason(reason)
  {
    exc_type="Exception";
  }

  
  void exception::PrintWhat() const
  {
    PrintWhat(std::cerr);
  }
  
  void exception::PrintWhat(std::ostream& where) const
  {
    where << what();
  }
  
  const char* exception::what() const throw()
  {
    std::stringstream ss;
    ss << exc_type << " due to " << m_reason << std::endl << "Please report code "
	  << m_code << " occured at "
       << m_file << ":" << m_line << "." << std::endl;
    return ss.str().c_str();
  }
  
  void exception::SetFile(const std::string& file)
  {
    m_file = file;
  }
    
  void exception::SetLine(const std::string& line)
  {
    m_line = line;
  }
  
  void exception::SetReason(const std::string& reason)
  {
    m_reason = reason;
  }
  
  const std::string& exception::GetReason() const
  {
    return m_reason;
  }
  
}
