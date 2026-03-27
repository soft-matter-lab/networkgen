; ============================================================================
; NetworkGen Installer Script
; Built with Inno Setup 6+  (https://jrsoftware.org/isinfo.php)
;
; HOW TO BUILD
; ------------
; 1. Install Inno Setup 6 from https://jrsoftware.org/isdl.php
; 2. Place this file in the root of the repository alongside the
;    NetworkGen/ source folder.
; 3. Open this file in the Inno Setup IDE and press Compile (Ctrl+F9),
;    OR run from the command line:
;       "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" NetworkGen.iss
; 4. The installer exe is written to installer_output\
;
; RELEASE WORKFLOW (GitHub Actions)
; -----------------------------------
; The CI workflow reads AppVersion from the VERSION file at repo root.
; To bump the version, edit VERSION (e.g. "1.3") before merging to main.
; The release job runs:
;   ISCC.exe /DAppVersion=%VERSION% NetworkGen.iss
; which overrides the default below.
; ============================================================================

; ── Version ──────────────────────────────────────────────────────────────────
; Override on the command line with /DAppVersion=X.Y
#ifndef AppVersion
  #define AppVersion "0.1"
#endif

; ── App metadata ─────────────────────────────────────────────────────────────
#define AppName        "NetworkGen"
#define AppPublisher   "Soft Matter Lab"
#define AppURL         "https://soft-matter-lab.github.io/networkgen"
#define AppExeName     ""          ; no exe — MATLAB tool, no launcher needed
#define SourceDir      "NetworkGen"   ; folder containing the .m files

[Setup]
AppId={{B3A7C2D1-4F8E-4A1B-9C3D-2E5F6A7B8C9D}
AppName={#AppName}
AppVersion={#AppVersion}
AppPublisher={#AppPublisher}
AppPublisherURL={#AppURL}
AppSupportURL={#AppURL}
AppUpdatesURL={#AppURL}/releases

; Default install location: Documents\MATLAB\NetworkGen
; This matches MATLAB's default userpath on Windows.
DefaultDirName={userdocs}\MATLAB\NetworkGen
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes

; Output
OutputDir=installer_output
OutputBaseFilename=NetworkGen_v{#AppVersion}_Setup

; Appearance
WizardStyle=modern

; Compression
Compression=lzma2/ultra64
SolidCompression=yes
InternalCompressLevel=ultra64

; Minimum Windows version (Windows 10)
MinVersion=10.0.17763

; Do not require admin rights — installs per-user into Documents
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

; ── Files ────────────────────────────────────────────────────────────────────
[Files]
; Copy all .m files (and subfolders) from the NetworkGen source folder
Source: "{#SourceDir}\*"; \
    DestDir: "{app}"; \
    Flags: ignoreversion recursesubdirs createallsubdirs

; ── Icons ────────────────────────────────────────────────────────────────────
; No desktop shortcut or Start Menu entry — this is a MATLAB toolbox,
; not a standalone application.  Users access it through MATLAB.

; ── Registry (optional) ──────────────────────────────────────────────────────
; Record the install path so other tools can find NetworkGen
[Registry]
Root: HKCU; \
    Subkey: "Software\SoftMatterLab\NetworkGen"; \
    ValueType: string; ValueName: "InstallPath"; \
    ValueData: "{app}"; \
    Flags: uninsdeletekey

; ── Pascal script ────────────────────────────────────────────────────────────
[Code]

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function GetMatlabStartupFile(): String;
// Returns the path to the user's MATLAB startup.m.
// MATLAB's default userpath on Windows is Documents\MATLAB.
var
  MatlabDir: String;
begin
  MatlabDir := ExpandConstant('{userdocs}\MATLAB');
  ForceDirectories(MatlabDir);
  Result := MatlabDir + '\startup.m';
end;


function StartupContainsMarker(StartupFile: String; Marker: String): Boolean;
// Returns True if startup.m already contains our marker string.
var
  Lines: TArrayOfString;
  i: Integer;
begin
  Result := False;
  if not FileExists(StartupFile) then Exit;
  if LoadStringsFromFile(StartupFile, Lines) then
    for i := 0 to GetArrayLength(Lines) - 1 do
      if Pos(Marker, Lines[i]) > 0 then
      begin
        Result := True;
        Break;
      end;
end;


procedure InjectStartupPath();
// Appends one addpath line to startup.m so MATLAB always finds NetworkGen.
// Safe to call multiple times — checks for the marker first.
var
  StartupFile : String;
  AppPath     : String;
  Marker      : String;
  Line        : String;
  F           : Integer;
begin
  StartupFile := GetMatlabStartupFile();
  AppPath     := ExpandConstant('{app}');
  Marker      := 'NetworkGen installer';
  Line        := 'addpath(genpath(''' + AppPath + ''')); % Added by ' + Marker;

  if StartupContainsMarker(StartupFile, Marker) then
    Exit;  // already present from a previous install

  // Open for append (create if missing)
  F := FileOpen(StartupFile, fmOpenWrite);
  if F = -1 then
    F := FileCreate(StartupFile);

  FileSeek(F, 0, 2);    // seek to end of file
  FileWrite(F, #13#10 + Line + #13#10);
  FileClose(F);
end;


procedure RemoveStartupPath();
// Removes the NetworkGen line from startup.m during uninstall.
var
  StartupFile : String;
  Lines       : TArrayOfString;
  NewLines    : TArrayOfString;
  i, j        : Integer;
  Marker      : String;
begin
  StartupFile := GetMatlabStartupFile();
  Marker      := 'NetworkGen installer';

  if not FileExists(StartupFile) then Exit;
  if not LoadStringsFromFile(StartupFile, Lines) then Exit;

  j := 0;
  SetArrayLength(NewLines, GetArrayLength(Lines));

  for i := 0 to GetArrayLength(Lines) - 1 do
  begin
    if Pos(Marker, Lines[i]) = 0 then  // keep lines that don't contain marker
    begin
      NewLines[j] := Lines[i];
      j := j + 1;
    end;
  end;

  SetArrayLength(NewLines, j);
  SaveStringsToFile(StartupFile, NewLines, False);
end;


// ---------------------------------------------------------------------------
// Event hooks
// ---------------------------------------------------------------------------

procedure CurStepChanged(CurStep: TSetupStep);
begin
  if CurStep = ssPostInstall then
    InjectStartupPath();
end;


procedure CurUninstallStepChanged(CurUninstallStep: TUninstallStep);
begin
  if CurUninstallStep = usPostUninstall then
    RemoveStartupPath();
end;


// ---------------------------------------------------------------------------
// Upgrade detection — offer to uninstall old version first
// ---------------------------------------------------------------------------

function InitializeSetup(): Boolean;
var
  OldPath     : String;
  Uninstaller : String;
  ResultCode  : Integer;
  Msg         : String;
begin
  Result := True;

  // Check registry for an existing install
  if RegQueryStringValue(HKCU,
      'Software\SoftMatterLab\NetworkGen',
      'InstallPath', OldPath) then
  begin
    if DirExists(OldPath) then
    begin
      Msg := 'An existing NetworkGen installation was found at:' + #13#10
           + OldPath + #13#10 + #13#10
           + 'It is recommended to remove it before installing the new version.' + #13#10 + #13#10
           + 'Remove it now?';

      if MsgBox(Msg, mbConfirmation, MB_YESNO) = IDYES then
      begin
        Uninstaller := OldPath + '\unins000.exe';
        if FileExists(Uninstaller) then
        begin
          Exec(Uninstaller, '/SILENT', '', SW_SHOW, ewWaitUntilTerminated, ResultCode);
        end else
        begin
          // No uninstaller found — remove directory and clean startup.m
          DelTree(OldPath, True, True, True);
          RemoveStartupPath();
        end;
      end;
    end;
  end;
end;
