if {[catch {package require Tcl 8.5}]} return
set script ""
if {![info exists ::env(TREECTRL_LIBRARY)]
    && [file exists [file join $dir treectrl.tcl]]} {
    append script "[list set ::treectrl_library $dir]\n"
}
append script [list load [file join $dir treectrl24.dll] treectrl]
package ifneeded treectrl 2.4.1 $script
