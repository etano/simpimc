#include "WriteClass.h"

void Write::DoEvent()
{
  vector<Event*>::iterator iter;
  for (iter=events.begin(); iter!=events.end(); ++iter) {
    (*iter)->Write();
  }
}
