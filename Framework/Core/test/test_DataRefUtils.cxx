// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <TObject.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMessage.h>
#include "Framework/RootSerializationSupport.h"
#include "Framework/DataRefUtils.h"
#include <catch_amalgamated.hpp>

#include <memory>

using namespace o2::framework;

TEST_CASE("PureRootTest")
{
  TBufferFile buffer(TBuffer::kWrite);
  TObjString s("test");
  buffer.WriteObject(&s);

  TBufferFile buffer2(TBuffer::kRead, buffer.BufferSize(), buffer.Buffer(), false);
  buffer2.SetReadMode();
  buffer2.InitMap();
  TClass* storedClass = buffer2.ReadClass();
  // ReadClass advances the buffer, so we need to reset it.
  buffer2.SetBufferOffset(0);
  buffer2.ResetMap();
  REQUIRE(storedClass != nullptr);
  auto* outS = (TObjString*)buffer2.ReadObjectAny(storedClass);
  REQUIRE(outS != nullptr);
  REQUIRE(outS->GetString() == "test");
}

// Simple test to do root deserialization.
TEST_CASE("TestRootSerialization")
{
  DataRef ref;
  auto transport = fair::mq::TransportFactory::CreateTransportFactory("zeromq");
  auto msg = transport->CreateMessage(4096);
  FairOutputTBuffer tm(*msg);
  auto sOrig = std::make_unique<TObjString>("test");
  tm << sOrig.get();
  o2::header::DataHeader dh;
  dh.payloadSerializationMethod = o2::header::gSerializationMethodROOT;
  ref.payload = (char*)msg->GetData();
  dh.payloadSize = (size_t)msg->GetSize();
  ref.header = reinterpret_cast<char const*>(&dh);

  // Check by using the same type
  auto s = DataRefUtils::as<TObjString>(ref);
  REQUIRE(s.get() != nullptr);
  REQUIRE(std::string(s->GetString().Data()) == "test");
  REQUIRE(std::string(s->GetName()) == "test");

  // Check by using the base type.
  auto o = DataRefUtils::as<TObject>(ref);
  REQUIRE(o.get() != nullptr);
  REQUIRE(std::string(o->GetName()) == "test");
}

// Simple test for ROOT container deserialization.
TEST_CASE("TestRootContainerSerialization")
{
  DataRef ref;
  TMessage* tm = new TMessage(kMESS_OBJECT);
  TObjArray container;
  // the original container is explicitly non-owning, its owning
  // flag is preserved during serialization
  container.SetOwner(false);
  auto object = std::make_unique<TObjString>("test");
  container.Add(object.get());
  tm->WriteObject(&container);
  o2::header::DataHeader dh;
  dh.payloadSerializationMethod = o2::header::gSerializationMethodROOT;
  ref.payload = tm->Buffer();
  dh.payloadSize = tm->BufferSize();
  ref.header = reinterpret_cast<char const*>(&dh);

  auto s = DataRefUtils::as<TObjArray>(ref);
  REQUIRE(s.get() != nullptr);
  REQUIRE(s->GetEntries() == 1);
  REQUIRE(std::string(s->At(0)->GetName()) == "test");
  // the extracted object must be owning to avoid memory leaks
  // the get method takes care of this
  REQUIRE(s->IsOwner());
}
