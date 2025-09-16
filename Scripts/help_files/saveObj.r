function (FileName = "", ObjName = "", DataDir = file.path(WRKDIR, 
    "Data")) 
{
    assign(FileName, get(ObjName))
    save(list = FileName, file = file.path(DataDir, FileName))
    rm(list = FileName)
}
