#include "StorkApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<StorkApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

StorkApp::StorkApp(InputParameters parameters) : MooseApp(parameters)
{
  StorkApp::registerAll(_factory, _action_factory, _syntax);
}

StorkApp::~StorkApp() {}

void
StorkApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"StorkApp"});
  Registry::registerActionsTo(af, {"StorkApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
StorkApp::registerApps()
{
  registerApp(StorkApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
StorkApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  StorkApp::registerAll(f, af, s);
}
extern "C" void
StorkApp__registerApps()
{
  StorkApp::registerApps();
}
