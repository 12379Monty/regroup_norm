function (FileName = "", ObjName = "", DataDir = file.path(WRKDIR, 
    "Data")) 
{
    load(file.path(DataDir, FileName))
    assign(ObjName, get(FileName), pos = 1)
    rm(list = FileName)
}
