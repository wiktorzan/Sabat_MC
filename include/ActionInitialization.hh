/// \file ActionInitialization.hh
/// \brief Definition of the ActionInitialization class (Mandatory)

#ifndef ACTION_INITIALIZATION_HH
#define ACTION_INITIALIZATION_HH

#include <G4VUserActionInitialization.hh>
#include <string>

/// Action initialization class.
///
class ActionInitialization : public G4VUserActionInitialization
{
public:
  ///Constructor
  ActionInitialization();
  ///Destructor
  ~ActionInitialization();

  ///Basic interface intherited from the base class G4VUserActionInitialization

  ///Register User Actions for the workers
  void Build() const override;

  ///Registers User Actions for the master (only Run Action)
  void BuildForMaster() const override;

  void SetTimeAndSeedAdd(std::string add) {timeAndSeedAdd = add;};
private:
  std::string timeAndSeedAdd = "";
};

#endif
