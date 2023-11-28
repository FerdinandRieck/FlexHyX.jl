#-- Datenstruktur
abstract type flexhyx end

abstract type Knoten <: flexhyx end
abstract type Kante <: flexhyx end

abstract type Strom_Knoten <: Knoten end
abstract type Gas_Knoten <: Knoten end
abstract type Temp_Knoten <: Knoten end
abstract type Wasser_Knoten <: Knoten end

abstract type Strom_Kante <: Kante end
abstract type Gas_Kante <: Kante end
abstract type Temp_Kante <: Kante end
abstract type Wasser_Kante <: Kante end
#----------------------------------------