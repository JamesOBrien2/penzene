; Inno Setup script: iscc /DVersion=x.y.z /DSource=dist\penzene cmake\penzene.iss
[Setup]
AppName=Penzene
AppVersion={#Version}
AppPublisher=Penzene contributors
AppPublisherURL=https://github.com/JamesOBrien2/penzene
DefaultDirName={autopf}\Penzene
DefaultGroupName=Penzene
UninstallDisplayIcon={app}\penzene.exe
OutputDir=.
OutputBaseFilename=penzene-windows-x64-setup
LicenseFile=..\LICENSE
Compression=lzma2
SolidCompression=yes
ArchitecturesInstallIn64BitMode=x64compatible
ArchitecturesAllowed=x64compatible
PrivilegesRequiredOverridesAllowed=dialog
ChangesAssociations=yes

[Files]
Source: "{#Source}\*"; DestDir: "{app}"; Flags: recursesubdirs ignoreversion

[Icons]
Name: "{group}\Penzene"; Filename: "{app}\penzene.exe"
Name: "{autodesktop}\Penzene"; Filename: "{app}\penzene.exe"; Tasks: desktopicon

[Tasks]
Name: desktopicon; Description: "Create a desktop shortcut"; Flags: unchecked

[Registry]
Root: HKA; Subkey: "Software\Classes\.penz"; ValueType: string; ValueData: "Penzene.Document"; Flags: uninsdeletevalue
Root: HKA; Subkey: "Software\Classes\Penzene.Document"; ValueType: string; ValueData: "Penzene document"; Flags: uninsdeletekey
Root: HKA; Subkey: "Software\Classes\Penzene.Document\DefaultIcon"; ValueType: string; ValueData: "{app}\penzene.exe,0"
Root: HKA; Subkey: "Software\Classes\Penzene.Document\shell\open\command"; ValueType: string; ValueData: """{app}\penzene.exe"" ""%1"""

; Penzene Drawing Object: Word and PowerPoint embed it, and double-click edits it here (#229, src/OleServer.cpp).
; OLE's default handler (ole32.dll) draws the cached picture; COM starts penzene.exe -Embedding to edit.
Root: HKA; Subkey: "Software\Classes\Penzene.Drawing"; ValueType: string; ValueData: "Penzene Drawing"; Flags: uninsdeletekey
Root: HKA; Subkey: "Software\Classes\Penzene.Drawing\CLSID"; ValueType: string; ValueData: "{{F6858801-14D9-488C-B79A-8C7D07452862}"
Root: HKA; Subkey: "Software\Classes\Penzene.Drawing\Insertable"; ValueType: none
Root: HKA; Subkey: "Software\Classes\CLSID\{{F6858801-14D9-488C-B79A-8C7D07452862}"; ValueType: string; ValueData: "Penzene Drawing"; Flags: uninsdeletekey
Root: HKA; Subkey: "Software\Classes\CLSID\{{F6858801-14D9-488C-B79A-8C7D07452862}\ProgID"; ValueType: string; ValueData: "Penzene.Drawing"
Root: HKA; Subkey: "Software\Classes\CLSID\{{F6858801-14D9-488C-B79A-8C7D07452862}\LocalServer32"; ValueType: string; ValueData: """{app}\penzene.exe"""
Root: HKA; Subkey: "Software\Classes\CLSID\{{F6858801-14D9-488C-B79A-8C7D07452862}\InprocHandler32"; ValueType: string; ValueData: "ole32.dll"
Root: HKA; Subkey: "Software\Classes\CLSID\{{F6858801-14D9-488C-B79A-8C7D07452862}\DefaultIcon"; ValueType: string; ValueData: "{app}\penzene.exe,0"
Root: HKA; Subkey: "Software\Classes\CLSID\{{F6858801-14D9-488C-B79A-8C7D07452862}\Insertable"; ValueType: none
Root: HKA; Subkey: "Software\Classes\CLSID\{{F6858801-14D9-488C-B79A-8C7D07452862}\MiscStatus"; ValueType: string; ValueData: "0"
Root: HKA; Subkey: "Software\Classes\CLSID\{{F6858801-14D9-488C-B79A-8C7D07452862}\AuxUserType\2"; ValueType: string; ValueData: "Penzene"
Root: HKA; Subkey: "Software\Classes\CLSID\{{F6858801-14D9-488C-B79A-8C7D07452862}\AuxUserType\3"; ValueType: string; ValueData: "Penzene"
Root: HKA; Subkey: "Software\Classes\CLSID\{{F6858801-14D9-488C-B79A-8C7D07452862}\Verb\0"; ValueType: string; ValueData: "&Edit,0,2"
Root: HKA; Subkey: "Software\Classes\CLSID\{{F6858801-14D9-488C-B79A-8C7D07452862}\Verb\1"; ValueType: string; ValueData: "&Open,0,2"
; The picture it offers: CF_ENHMETAFILE and CF_METAFILEPICT, content aspect, get only.
Root: HKA; Subkey: "Software\Classes\CLSID\{{F6858801-14D9-488C-B79A-8C7D07452862}\DataFormats\GetSet\0"; ValueType: string; ValueData: "14,1,64,1"
Root: HKA; Subkey: "Software\Classes\CLSID\{{F6858801-14D9-488C-B79A-8C7D07452862}\DataFormats\GetSet\1"; ValueType: string; ValueData: "3,1,32,1"

[Run]
Filename: "{app}\penzene.exe"; Description: "Launch Penzene"; Flags: nowait postinstall skipifsilent
