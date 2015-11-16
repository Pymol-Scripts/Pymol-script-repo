if {[catch {package require Tcl 8.5}]} return
set script "load \"[file join $dir shellicon24.dll]\" shellicon"
package ifneeded shellicon 2.4.1 $script
